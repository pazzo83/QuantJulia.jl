using StatsFuns

type FdmLinearOpLayout{I <: Integer}
  size::I
  dim::Vector{I}
  spacing::Vector{I}
end

function FdmLinearOpLayout{I <: Integer}(dim::Vector{I})
  spacing = ones(Int, length(dim))
  spacing[2:end] = cumprod(dim[1:end-1])
  sz = spacing[end] * dim[end]

  return FdmLinearOpLayout(sz, dim, spacing)
end

function get_layout_from_meshers{F1D <: Fdm1DMesher}(mesherVec::Vector{F1D})
  dim = zeros(Int, length(mesherVec))
  for i = 1:length(dim)
    dim[i] = mesherVec[i].size
  end

  return FdmLinearOpLayout(dim)
end

type FdmBoundaryConditionSet end

type FdmDividendHandler{F <: FdmMesher, I <: Integer} <: StepCondition
  x::Vector{Float64}
  dividendTimes::Vector{Float64}
  dividendDates::Vector{Date}
  dividends::Vector{Float64}
  mesher::F
  equityDirection::I
end

function FdmDividendHandler(schedule::DividendSchedule, mesher::FdmMesher, refDate::Date, dc::DayCount, equityDirection::Int)
  schedLength = length(schedule.dividends)
  dividends = zeros(schedLength)
  dividendDates = Vector{Date}(schedLength)
  dividendTimes = zeros(schedLength)
  x = zeros(mesher.layout.dim[equityDirection])

  for i = 1:schedLength
    dividends = amount(schedule.dividends[i])
    dividendDates = date(schedule.dividends[i])
    dividendTimes = year_fraction(dc, refDate, date(schedule.dividends[i]))
  end

  tmp = get_locations(mesher, equityDirection)
  spacing = mesher.layout.spacing[equityDirection]

  for i = 1:length(x)
    iter = ((i - 1) * spacing) + 1
    x[i] = exp(tmp[iter])
  end

  return FdmDividendHandler(x, dividendTimes, dividendDates, dividends, mesher, equityDirection)
end

type FdmAmericanStepCondition{F <: FdmMesher, C <: FdmInnerValueCalculator} <: StepCondition
  mesher::F
  calculator::C
end

type FdmBermudanStepCondition{F <: FdmMesher, C <: FdmInnerValueCalculator} <: StepCondition
  mesher::F
  calculator::C
  exerciseTimes::Vector{Float64}
end

function FdmBermudanStepCondition(exerciseDates::Vector{Date}, refDate::Date, dc::DayCount, mesher::FdmMesher, calculator::FdmInnerValueCalculator)
  exerciseTimes = zeros(length(exerciseDates))
  for i = 1:length(exerciseTimes)
    exerciseTimes[i] = year_fraction(dc, refDate, exerciseDates[i])
  end

  return FdmBermudanStepCondition(mesher, calculator, exerciseTimes)
end

type FdmStepConditionComposite{C <: StepCondition}
  stoppingTimes::Vector{Float64}
  conditions::Vector{C}
end

# Constructors
function vanilla_FdmStepConditionComposite(cashFlow::DividendSchedule, exercise::Exercise, mesher::FdmMesher, calculator::FdmInnerValueCalculator,
                                          refDate::Date, dc::DayCount)
  stoppingTimes = Vector{Float64}()
  stepConditions = Vector{StepCondition}()

  if length(cashFlow.dividends) > 0
    dividendCondition = FdmDividendHandler(cashFlow, mesher, refDate, dc, 1)
    push!(stepConditions, dividendCondition)
    append!(stoppingTimes, dividendCondition.dividendTimes)
  end

  if isa(exercise, AmericanExercise)
    push!(stepConditions, FdmAmericanStepCondition(mesher, calculator))
  elseif isa(exercise, BermudanExercise)
    bermudanCondition = FdmBermudanStepCondition(exercise.dates, refDate, dc, mesher, calculator)
    push!(stepConditions, bermudanCondition)
    append!(stoppingTimes, bermudanCondition.exerciseTimes)
  end

  return FdmStepConditionComposite(stoppingTimes, stepConditions)
end

type FdmMesherComposite{FM1D <: Fdm1DMesher} <: FdmMesher
  layout::FdmLinearOpLayout
  meshers::Vector{FM1D} # this could change
end

function FdmMesherComposite{F1D <: Fdm1DMesher}(xmesher::F1D, ymesher::F1D)
  meshers = F1D[xmesher, ymesher]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite(layout, meshers)
end

function iter_coords!(coord::Vector{Int}, dims::Vector{Int})
  for i = 1:length(dims)
    coord[i] += 1
    if coord[i] == dim[i] + 1
      coord[i] = 1
    else
      break
    end
  end

  return coord
end

function get_locations(mesher::FdmMesherComposite, direction::Int)
  coords = ones(length(mesher.layout.dim))
  retVal = zeros(mesher.layout.size)
  for i = 1:length(retVal)
    retVal[i] = mesher.meshers[i].locations[coord[direction]]
    iter_coords!(coords, mesher.layout.dim)
  end

  return retVal
end

type FdmAffineModelSwapInnerValue{T1 <: TermStructure, T2 <: TermStructure, M1 <: Model, M2 <: Model, FM <: FdmMesher, I <: Integer} <: FdmInnerValueCalculator
  disTs::T1
  fwdTs::T2
  disModel::M1
  fwdModel::M2
  swap::VanillaSwap
  exerciseDates::Dict{Float64, Date}
  mesher::FM
  direction::I
end

function FdmAffineModelSwapInnerValue(disModel::Model, fwdModel::Model, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FdmMesher, direction::Int)
  newSwap = VanillaSwap(swap.swapT, swap.nominal, swap.fixedSchedule, swap.fixedRate, swap.fixedDayCount, swap.iborIndex, swap.spread, swap.floatSchedule,
                        swap.floatDayCount, swap.pricingEngine, swap.paymentConvention)

  return FdmAffineModelSwapInnerValue(disModel.ts, fwdModel.ts, disModel, fwdModel, newSwap, exerciseDates, mesher, direction)
end

type Hundsdorfer <: FdmSchemeDescType end

type FdmSchemeDesc{F <: FdmSchemeDescType}
  schemeType::F
  theta::Float64
  mu::Float64
end

FdmSchemeDesc(t::Hundsdorfer) = FdmSchemeDesc(t, 0.5 + sqrt(3.0) / 6.0, 0.5)

## Meshers ##
type FdmSimpleProcess1dMesher{I <: Integer, P <: StochasticProcess1D} <: Fdm1DMesher
  size::I
  process::P
  maturity::Float64
  tAvgSteps::I
  epsilon::Float64
  mandatoryPoint::Float64
  locations::Vector{Float64}
  dplus::Vector{Float64}
  dminus::Vector{Float64}
end

function FdmSimpleProcess1dMesher(sz::Int, process::StochasticProcess1D, maturity::Float64, tAvgSteps::Int, _eps::Float64, mandatoryPoint::Float64 = -1.0)
  locations = zeros(sz)
  dminus = zeros(sz)
  dplus = zeros(sz)
  mp = mandatoryPoint == -1.0 ? process.x0 : mandatoryPoint
  for l = 1:tAvgSteps
    t = (maturity * l) / tAvgSteps

    qMin = min(process.x0, evolve(process, 0.0, process.x0, t, norminvcdf(_eps)))
    qMax = max(max(process.x0, mp), evolve(process, 0.0, process.x0, t, norminvcdf(1.0 - _eps)))

    dp = (1.0 - 2.0 * _eps) / (sz - 1)
    p = _eps
    locations[1] += qMin

    for i = 2:sz - 1
      p += dp
      locations[i] += evolve(process, 0.0, process.x0, t, norminvcdf(p))
    end

    locations[end] += qMax
  end

  locations /= tAvgSteps
  for i = 1:sz - 1
    dminus[i + 1] = dplus[i] = locations[i + 1] - locations[i]
  end

  dplus[end] = dminus[1] = -1.0

  return FdmSimpleProcess1dMesher(sz, process, maturity, tAvgSteps, _eps, mandatoryPoint, locations, dplus, dminus)
end

type FdmSolverDesc{F <: FdmMesher, C <: FdmInnerValueCalculator, I <: Integer}
  mesher::F
  bcSet::FdmBoundaryConditionSet
  condition::FdmStepConditionComposite
  calculator::C
  maturity::Float64
  timeSteps::I
  dampingSteps::I
end

type FdmG2Solver <: LazyObject
  lazyMixin::LazyMixin
  model::G2
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc

  FdmG2Solver(model::G2, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc) = new(LazyMixin(), model, solverDesc, schemeDesc)
end
