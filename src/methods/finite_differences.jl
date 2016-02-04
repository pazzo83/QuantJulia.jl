using QuantJulia.Math
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

function neighborhood{I <: Integer}(mesherLayout::FdmLinearOpLayout, idx::I, coords::Vector{I}, i::I, offset::I)
  myIndex = idx - (coords[i] - 1) * mesherLayout.spacing[i]
  coorOffset = (coords[i] - 1) + offset

  if coorOffset < 0
    coorOffset = -coorOffset
  elseif coorOffset >= mesherLayout.dim[i]
    coorOffset = 2 * (mesherLayout.dim[i] - 1) - coorOffset
  end

  return myIndex + coorOffset * mesherLayout.spacing[i]
end

function neighborhood{I <: Integer}(mesherLayout::FdmLinearOpLayout, idx::I, coords::Vector{I}, i1::I, offset1::I, i2::I, offset2::I)
  myIndex = idx - (coords[i1] - 1) * mesherLayout.spacing[i1] - (coords[i2] - 1) * mesherLayout.spacing[i2]
  coorOffset1 = (coords[i1] - 1) + offset1

  if coorOffset1 < 0
    coorOffset1 = -coorOffset1
  elseif coorOffset1 >= mesherLayout.dim[i1]
    coorOffset1 = 2 * (mesherLayout.dim[i1] - 1) - coorOffset1
  end

  coorOffset2 = (coords[i2] - 1) + offset2
  if coorOffset2 < 0
    coorOffset2 = -coorOffset2
  elseif coorOffset2 >= mesherLayout.dim[i2]
    coorOffset2 = 2 * (mesherLayout.dim[i2] - 1) - coorOffset2
  end

  return myIndex + coorOffset1 * mesherLayout.spacing[i1] + coorOffset2 * mesherLayout.spacing[i2]
end

type FdmBoundaryConditionSet
  conditions::Vector{BoundaryCondition}
end

FdmBoundaryConditionSet() = FdmBoundaryConditionSet(Vector{BoundaryCondition}(0))

function set_time!(bcSet::FdmBoundaryConditionSet, t::Float64)
  for cond in bcSet.conditions
    set_time!(cond, t)
  end

  return bcSet
end

function apply_before_applying!(bcSet::FdmBoundaryConditionSet, op::FdmLinearOpComposite)
  for cond in bcSet.conditions
    apply_before_applying!(cond, t)
  end

  return bcSet
end

function apply_after_applying!(bcSet::FdmBoundaryConditionSet, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_after_applying!(cond, t)
  end

  return bcSet
end

function apply_before_solving!(bcSet::FdmBoundaryConditionSet, op::FdmLinearOpComposite, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_before_solving!(cond, op, a)
  end

  return bcSet
end

function apply_after_solving!(bcSet::FdmBoundaryConditionSet, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_after_solving!(cond, t)
  end

  return bcSet
end

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

function apply_to!(cond::FdmBermudanStepCondition, a::Vector{Float64}, t::Float64)
  if findfirst(cond.exerciseTimes, t) != 0
    layout = cond.mesher.layout
    dims = length(layout.dim)
    coords = ones(Int, dims)

    locations = zeros(dims)

    for i = 1:layout.size
      for j = 1:dims
        locations[dims] = get_location(cond.mesher, coords, j)
      end

      innerValue = inner_value(cond.calculator, coords, i, t)
      if innerValue > a[i]
        a[i] = innerValue
      end

      iter_coords!(coords, cond.mesher.layout.dim)
    end
  end

  return cond, a
end

type FdmSnapshotCondition <: StepCondition
  t::Float64
  a::Vector{Float64}
end

FdmSnapshotCondition(t::Float64) = FdmSnapshotCondition(t, Vector{Float64}(0))

function apply_to!(cond::FdmSnapshotCondition, a::Vector{Float64}, t::Float64)
  if cond.t == t
    cond.a = a
  end

  return cond
end

type FdmStepConditionComposite{C <: StepCondition} <: StepCondition
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

function join_conditions_FdmStepConditionComposite(c1::FdmSnapshotCondition, c2::FdmStepConditionComposite)
  stoppingTimes = c2.stoppingTimes
  push!(stoppingTimes, c1.t)

  conditions = Vector{StepCondition}(2)
  conditions[1] = c2
  conditions[2] = c1

  return FdmStepConditionComposite(stoppingTimes, conditions)
end

function apply_to!(cond::FdmStepConditionComposite, a::Vector{Float64}, t::Float64)
  for c in cond.conditions
    apply_to!(c, a, t)
  end

  return cond, a
end

type FdmMesherComposite{FM1D <: Fdm1DMesher} <: FdmMesher
  layout::FdmLinearOpLayout
  meshers::Vector{FM1D} # this could change
end

# Constructors
function FdmMesherComposite{F1D <: Fdm1DMesher}(mesh::F1D)
  meshers = F1D[mesh]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite(layout, meshers)
end

function FdmMesherComposite{F1D <: Fdm1DMesher}(xmesher::F1D, ymesher::F1D)
  meshers = F1D[xmesher, ymesher]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite(layout, meshers)
end

function iter_coords!(coord::Vector{Int}, dims::Vector{Int})
  for i = 1:length(dims)
    coord[i] += 1
    if coord[i] == dims[i] + 1
      coord[i] = 1
    else
      break
    end
  end

  return coord
end

function get_locations(mesher::FdmMesherComposite, direction::Int)
  coords = ones(Int, length(mesher.layout.dim))
  retVal = zeros(mesher.layout.size)
  for i = 1:length(retVal)
    retVal[i] = mesher.meshers[direction].locations[coords[direction]]
    iter_coords!(coords, mesher.layout.dim)
  end

  return retVal
end

get_location{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].locations[coords[direction]]

get_dminus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dminus[coords[direction]]
get_dplus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dplus[coords[direction]]

type FdmAffineModelTermStructure{I <: Integer, B <: BusinessCalendar, DC <: DayCount, A <: AffineModel} <: YieldTermStructure
  settlement_days::I
  referenceDate::Date
  calendar::B
  dc::DC
  modelReferenceDate::Date
  model::A
  r::Vector{Float64}
  t::Float64
end

FdmAffineModelTermStructure{B <: BusinessCalendar, DC <: DayCount, A <: AffineModel}(referenceDate::Date, cal::B, dc::DC, modelReferenceDate::Date, model::A, r::Vector{Float64}) =
                            FdmAffineModelTermStructure{Int, B, DC, A}(0, referenceDate, cal, dc, modelReferenceDate, model, r, year_fraction(dc, modelReferenceDate, referenceDate))

discount_impl(ts::FdmAffineModelTermStructure, T::Float64) = discount_bond(ts.model, ts.t, T + ts.t, ts.r)

set_variable(ts::FdmAffineModelTermStructure, r::Vector{Float64}) = ts.r = r

type FdmAffineModelSwapInnerValue{M1 <: Model, M2 <: Model, FM <: FdmMesher, I <: Integer} <: FdmInnerValueCalculator
  disModel::M1
  fwdModel::M2
  swap::VanillaSwap
  exerciseDates::Dict{Float64, Date}
  mesher::FM
  direction::I
  disTs::FdmAffineModelTermStructure
  fwdTs::FdmAffineModelTermStructure

  function FdmAffineModelSwapInnerValue{M1, M2, FM, I}(disModel::M1, fwdModel::M2, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FM, direction::I)
    idx = swap.iborIndex
    newIbor = IborIndex(idx.familyName, idx.tenor, idx.fixingDays, idx.currency, idx.fixingCalendar, idx.convention, idx.endOfMonth, idx.dc)
    newSwap = VanillaSwap(swap.swapT, swap.nominal, swap.fixedSchedule, swap.fixedRate, swap.fixedDayCount, newIbor, swap.spread, swap.floatSchedule,
                          swap.floatDayCount, swap.pricingEngine, swap.paymentConvention)

    return new{M1, M2, FM, I}(disModel, fwdModel, newSwap, exerciseDates, mesher, direction)
  end
end

FdmAffineModelSwapInnerValue{M1 <: Model, M2 <: Model, FM <: FdmMesher, I <: Integer}(disModel::M1, fwdModel::M2, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FM, direction::I) =
                            FdmAffineModelSwapInnerValue{M1, M2, FM, I}(disModel, fwdModel, swap, exerciseDates, mesher, direction)

function get_state{I <: Integer}(::G2, calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, ::Float64)
  retVal = Vector{Float64}(2)
  retVal[1] = get_location(calc.mesher, coords, calc.direction)
  retVal[2] = get_location(calc.mesher, coords, calc.direction + 1)

  return retVal
end

get_state{I <: Integer}(model::HullWhite, calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, t::Float64) = [short_rate(get_dynamics(model), t, get_location(calc.mesher, coords, calc.direction))]

function inner_value{I <: Integer}(calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, i::I, t::Float64)
  iterExerciseDate = get(calc.exerciseDates, t, Date())
  disRate = get_state(calc.disModel, calc, coords, t)
  fwdRate = get_state(calc.fwdModel, calc, coords, t)

  if !isdefined(calc, :disTs) || iterExerciseDate != reference_date(calc.disTs)
    disc = calc.disModel.ts
    calc.disTs = FdmAffineModelTermStructure(iterExerciseDate, disc.calendar, disc.dc, reference_date(disc), calc.disModel, disRate)
    fwd = calc.fwdModel.ts
    calc.fwdTs = FdmAffineModelTermStructure(iterExerciseDate, fwd.calendar, fwd.dc, reference_date(fwd), calc.fwdModel, fwdRate)
    # probably should put this in something more logical
    calc.swap.iborIndex.ts = calc.fwdTs
  else
    # do some updating
    set_variable(calc.disTs, disRate)
    set_variable(calc.fwdTs, fwdRate)
  end

  npv = 0.0
  for j = 1:2
    for i = 1:length(calc.swap.legs[j].coupons)
      cf = calc.swap.legs[j].coupons[i]
      if isa(cf, Coupon)
        npv += accrual_start_date(cf) >= iterExerciseDate ? amount(cf) * discount(calc.disTs, date(cf)) : 0.0
      end
    end
    if j == 1
      npv *= -1.0
    end
  end
  if isa(calc.swap.swapT, Receiver)
    npv *= -1.0
  end
  return max(0.0, npv)
end

avg_inner_value{I <: Integer}(calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, i::I, t::Float64) = inner_value(calc, coords, i, t)

type Hundsdorfer <: FdmSchemeDescType end
type Douglas <: FdmSchemeDescType end

type FdmSchemeDesc{F <: FdmSchemeDescType}
  schemeType::F
  theta::Float64
  mu::Float64
end

FdmSchemeDesc(t::Hundsdorfer) = FdmSchemeDesc(t, 0.5 + sqrt(3.0) / 6.0, 0.5)
FdmSchemeDesc(t::Douglas) = FdmSchemeDesc(t, 0.5, 0.0)

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

type TripleBandLinearOp{I <: Integer, FM <: FdmMesher}
  direction::I
  mesher::FM
  i0::Vector{I}
  i2::Vector{I}
  reverseIndex::Vector{I}
  lower::Vector{Float64}
  _diag::Vector{Float64}
  upper::Vector{Float64}
end

function TripleBandLinearOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function FirstDerivativeOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    ## FirstDerivativeOp part
    hm = get_dminus(mesher, coords, direction)
    hp = get_dplus(mesher, coords, direction)

    zetam1 = hm * (hm + hp)
    zeta0 = hm * hp
    zetap1 = hp * (hm + hp)

    if coords[direction] == 1
      # upwinding scheme
      lower[i] = 0.0
      upper[i] = 1.0 / hp
      _diag[i] = -(upper[i])
    elseif coords[direction] == mesher.layout.dim[direction]
      # downwinding scheme
      _diag[i] = 1.0 / hm
      lower[i] = -_diag[i]
      upper[i] = 0.0
    else
      lower[i] = -hp / zetam1
      _diag[i] = (hp - hm) / zeta0
      upper[i] = hm / zetap1
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function SecondDerivativeOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    ## SecondDerivativeOp part
    hm = get_dminus(mesher, coords, direction)
    hp = get_dplus(mesher, coords, direction)

    zetam1 = hm * (hm + hp)
    zeta0 = hm * hp
    zetap1 = hp * (hm + hp)

    co = coords[direction]

    if co == 1 || co == mesher.layout.dim[direction]
      lower[i] = upper[i] = _diag[i] = 0.0
    else
      lower[i] = 2.0 / zetam1
      _diag[i] = -2.0 / zeta0
      upper[i] = 2.0 / zetap1
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function mult!{T <: Number}(trpBandLinOp::TripleBandLinearOp, u::Vector{T})
  # probably should protect for array dims here
  trpBandLinOp.lower = trpBandLinOp.lower .* u
  trpBandLinOp._diag = trpBandLinOp._diag .* u
  trpBandLinOp.upper = trpBandLinOp.upper .* u

  return trpBandLinOp
end

function add!(trpBandLinOp::TripleBandLinearOp, m::TripleBandLinearOp)
  trpBandLinOp.lower += m.lower
  trpBandLinOp._diag += m._diag
  trpBandLinOp.upper += m.upper

  return trpBandLinOp
end

function axpyb!(trpBandLinOp::TripleBandLinearOp, a::Vector{Float64}, x::TripleBandLinearOp, y::TripleBandLinearOp, b::Vector{Float64})
  sz = trpBandLinOp.mesher.layout.size

  if length(a) == 0
    if length(b) == 0
      trpBandLinOp._diag[1:sz] = y._diag
      trpBandLinOp.lower[1:sz] = y.lower
      trpBandLinOp.upper[1:sz] = y.upper
    else
      addB = length(b) > 0 ? b[1:sz] : b[1]
      trpBandLinOp._diag[1:sz] = y._diag + addB
      trpBandLinOp.lower[1:sz] = y.lower
      trpBandLinOp.upper[1:sz] = y.upper
    end
  elseif length(b) == 0
    # this might be improved
    trpBandLinOp._diag[1:sz] = y._diag + (length(a) > 0 ? a .* x._diag : a[1] * x._diag)
    trpBandLinOp.lower[1:sz] = y.lower + (length(a) > 0 ? a .* x.lower : a[1] * x.lower)
    trpBandLinOp.upper[1:sz] = y.upper + (length(a) > 0 ? a .* x.upper : a[1] * x.upper)
  else
    addB = length(b) > 0 ? b[1:sz] : b[1]
    trpBandLinOp._diag[1:sz] = y._diag + (length(a) > 0 ? a .* x._diag : a[1] * x._diag) + addB
    trpBandLinOp.lower[1:sz] = y.lower + (length(a) > 0 ? a .* x.lower : a[1] * x.lower)
    trpBandLinOp.upper[1:sz] = y.upper + (length(a) > 0 ? a .* x.upper : a[1] * x.upper)
  end

  return trpBandLinOp
end

function apply(trpBandLinOp::TripleBandLinearOp, r::Vector{Float64})
  idx = trpBandLinOp.mesher.layout
  length(r) == idx.size || error("inconsistent length of r")

  retArray = Vector{Float64}(length(r))

  for i = 1:idx.size
    retArray[i] = r[trpBandLinOp.i0[i]] * trpBandLinOp.lower[i] + r[i] * trpBandLinOp._diag[i] + r[trpBandLinOp.i2[i]] * trpBandLinOp.upper[i]
  end

  return retArray
end

function solve_splitting(trpBandLinOp::TripleBandLinearOp, r::Vector{Float64}, a::Float64, b::Float64)
  layout = trpBandLinOp.mesher.layout
  length(r) == layout.size || error("inconsistent length of r")

  retArray = Vector{Float64}(length(r))
  tmp = Vector{Float64}(length(r))

  # solving a tridiagonal system (we could use Julia built in f'ns here)
  rim1 = trpBandLinOp.reverseIndex[1]
  bet = 1.0 / (a * trpBandLinOp._diag[rim1] + b)
  bet != 0.0 || error("division by zero")

  retArray[rim1] = r[rim1] * bet

  for j = 2:layout.size
    ri = trpBandLinOp.reverseIndex[j]
    tmp[j] = a * trpBandLinOp.upper[rim1] * bet

    bet = b + a * (trpBandLinOp._diag[ri] - tmp[j] * trpBandLinOp.lower[ri])
    bet != 0.0 || error("division by zero")

    bet = 1.0 / bet

    retArray[ri] = (r[ri] - a * trpBandLinOp.lower[ri] * retArray[rim1]) * bet
    rim1 = ri
  end

  for j = layout.size - 1:-1:2
    retArray[trpBandLinOp.reverseIndex[j]] -= tmp[j + 1] * retArray[trpBandLinOp.reverseIndex[j + 1]]
  end

  retArray[trpBandLinOp.reverseIndex[1]] -= tmp[2] * retArray[trpBandLinOp.reverseIndex[2]]

  return retArray
end

type SecondOrderMixedDerivativeOp{I <: Integer, FD <: FdmMesher} <: NinePointLinearOp
  d1::I
  d2::I
  i00::Vector{I}
  i10::Vector{I}
  i20::Vector{I}
  i01::Vector{I}
  i21::Vector{I}
  i02::Vector{I}
  i12::Vector{I}
  i22::Vector{I}
  a00::Vector{Float64}
  a10::Vector{Float64}
  a20::Vector{Float64}
  a01::Vector{Float64}
  a11::Vector{Float64}
  a21::Vector{Float64}
  a02::Vector{Float64}
  a12::Vector{Float64}
  a22::Vector{Float64}
  mesher::FD
end

function SecondOrderMixedDerivativeOp(d1::Int, d2::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i00 = Vector{Int}(sz)
  i10 = Vector{Int}(sz)
  i20 = Vector{Int}(sz)
  i01 = Vector{Int}(sz)
  i21 = Vector{Int}(sz)
  i02 = Vector{Int}(sz)
  i12 = Vector{Int}(sz)
  i22 = Vector{Int}(sz)
  a00 = Vector{Float64}(sz)
  a10 = Vector{Float64}(sz)
  a20 = Vector{Float64}(sz)
  a01 = Vector{Float64}(sz)
  a11 = Vector{Float64}(sz)
  a21 = Vector{Float64}(sz)
  a02 = Vector{Float64}(sz)
  a12 = Vector{Float64}(sz)
  a22 = Vector{Float64}(sz)

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    # NinePointLinearOp part
    i10[i] = neighborhood(mesher.layout, i, coords, d2, -1)
    i01[i] = neighborhood(mesher.layout, i, coords, d1, -1)
    i21[i] = neighborhood(mesher.layout, i, coords, d1,  1)
    i12[i] = neighborhood(mesher.layout, i, coords, d2,  1)
    i00[i] = neighborhood(mesher.layout, i, coords, d1, -1, d2, -1)
    i20[i] = neighborhood(mesher.layout, i, coords, d1,  1, d2, -1)
    i02[i] = neighborhood(mesher.layout, i, coords, d1, -1, d2,  1)
    i22[i] = neighborhood(mesher.layout, i, coords, d1,  1, d2,  1)

    # SecondOrderMixedDerivativeOp part
    hm_d1 = get_dminus(mesher, coords, d1)
    hp_d1 = get_dplus(mesher, coords, d1)
    hm_d2 = get_dminus(mesher, coords, d2)
    hp_d2 = get_dplus(mesher, coords, d2)

    zetam1 = hm_d1 * (hm_d1 + hp_d1)
    zeta0 = hm_d1 * hp_d1
    zetap1 = hp_d1 * (hm_d1 + hp_d1)

    phim1 = hm_d2 * (hm_d2 + hp_d2)
    phi0 = hm_d2 * hp_d2
    phip1 = hp_d2 * (hm_d2 + hp_d2)

    c1 = coords[d1]
    c2 = coords[d2]

    if c1 == 1 && c2 == 1
      # lower left corner
      a00[i] = a01[i] = a02[i] = a10[i] = a20[i] = 0.0
      a11[i] = a22[i] = 1.0 / (hp_d1 * hp_d2)
      a21[i] = a12[i] = -a11[i]
    elseif c1 == mesher.layout.dim[d1] && c2 == 1
      # upper left corner
      a22[i] = a21[i] = a20[i] = a10[i] = a00[i] = 0.0
      a01[i] = a12[i] = 1.0 / (hm_d1 * hp_d2)
      a11[i] = a02[i] = -a01[i]
    elseif c1 == 1 && c2 == mesher.layout.dim[d2]
      # lower right corner
      a00[i] = a01[i] = a02[i] = a12[i] = a22[i] = 0.0
      a10[i] = a21[i] = 1.0 / (hp_d1 * hm_d2)
      a20[i] = a11[i] = -a10[i]
    elseif c1 == mesher.layout.dim[d1] && c2 == mesher.layout.dim[d2]
      # upper right corner
      a20[i] = a21[i] = a22[i] = a12[i] = a02[i] = 0.0
      a00[i] = a11[i] = 1.0 / (hm_d1 * hm_d2)
      a10[i] = a01[i] = -a00[i]
    elseif c1 == 1
      # lower side
      a00[i] = a01[i] = a02[i] = 0.0
      a10[i] = hp_d2 / (hp_d1 * phim1)
      a20[i] = -a10[i]
      a21[i] = (hp_d2 - hm_d2) / (hp_d1 * phi0)
      a11[i] = -a21[i]
      a22[i] = hm_d2 / (hp_d1 * phip1)
      a12[i] = -a22[i]
    elseif c1 == mesher.layout.dim[d1]
      # upper side
      a20[i] = a21[i] = a22[i] = 0.0
      a00[i] = hp_d2 / (hm_d1 * phim1)
      a10[i] = -a00[i]
      a11[i] = (hp_d2 - hm_d2) / (hm_d1 * phi0)
      a01[i] = -a11[i]
      a12[i] = hm_d2 / (hm_d1 * phip1)
      a02[i] = -a12[i]
    elseif c2 == 1
      # left side
      a00[i] = a10[i] = a20[i] = 0.0
      a01[i] = hp_d1 / (zetam1 * hp_d2)
      a02[i] = -a01[i]
      a12[i] = (hp_d1 - hm_d1) / (zeta0 * hp_d2)
      a11[i] = -a12[i]
      a22[i] = hm_d1 / (zetap1 * hp_d2)
      a21[i] = -a22[i]
    elseif c2 == mesher.layout.dim[d2]
      # right side
      a22[i] = a12[i] = a02[i] = 0.0
      a00[i] = hp_d1 / (zetam1 * hm_d2)
      a01[i] = -a00[i]
      a11[i] = (hp_d1 - hm_d1) / (zeta0 * hm_d2)
      a10[i] = -a11[i]
      a21[i] = hm_d1 / (zetap1 * hm_d2)
      a20[i] = -a21[i]
    else
      a00[i] = hp_d1 * hp_d2 / (zetam1 * phim1)
      a10[i] = -(hp_d1 - hm_d1) * hp_d2 / (zeta0 * phim1)
      a20[i] = -hm_d1 * hp_d2 / (zetap1 * phim1)
      a01[i] = -hp_d1 * (hp_d2 - hm_d2) / (zetam1 * phi0)
      a11[i] = (hp_d1 - hm_d1) * (hp_d2 - hm_d2) / (zeta0 * phi0)
      a21[i] = hm_d1 * (hp_d2 - hm_d2) / (zetap1 * phi0)
      a02[i] = -hp_d1 * hm_d2 / (zetam1 * phip1)
      a12[i] = hm_d2 * (hp_d1 - hm_d1) / (zeta0 * phip1)
      a22[i] = hm_d1 * hm_d2 / (zetap1 * phip1)
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return SecondOrderMixedDerivativeOp(d1, d2, i00, i10, i20, i01, i21, i02, i12, i22, a00, a10, a20, a01, a11, a21, a02, a12, a22, mesher)
end

function mult!{T <: Number}(ninePointLin::NinePointLinearOp, u::Vector{T})
  ninePointLin.a11 = ninePointLin.a11 .* u
  ninePointLin.a01 = ninePointLin.a01 .* u
  ninePointLin.a10 = ninePointLin.a10 .* u
  ninePointLin.a21 = ninePointLin.a21 .* u
  ninePointLin.a22 = ninePointLin.a22 .* u
  ninePointLin.a00 = ninePointLin.a00 .* u
  ninePointLin.a02 = ninePointLin.a02 .* u
  ninePointLin.a20 = ninePointLin.a20 .* u
  ninePointLin.a12 = ninePointLin.a12 .* u

  return ninePointLin
end

function apply(ninePointLin::NinePointLinearOp, u::Vector{Float64})
  length(u) == ninePointLin.mesher.layout.size || error("inconsistent length of r")

  retVal = ninePointLin.a00 .* u[ninePointLin.i00] + ninePointLin.a01 .* u[ninePointLin.i01] + ninePointLin.a02 .* u[ninePointLin.i02] + ninePointLin.a10 .* u[ninePointLin.i10] +
          ninePointLin.a11 .* u + ninePointLin.a12 .* u[ninePointLin.i12] + ninePointLin.a20 .* u[ninePointLin.i20] + ninePointLin.a21 .* u[ninePointLin.i21] +
          ninePointLin.a22 .* u[ninePointLin.i22]

  return retVal
end

type FdmG2Op{I <: Integer} <: FdmLinearOpComposite
  direction1::I
  direction2::I
  x::Vector{Float64}
  y::Vector{Float64}
  dxMap::TripleBandLinearOp
  dyMap::TripleBandLinearOp
  corrMap::SecondOrderMixedDerivativeOp
  mapX::TripleBandLinearOp
  mapY::TripleBandLinearOp
  model::G2
end

function FdmG2Op(mesher::FdmMesher, model::G2, direction1::Int, direction2::Int)
  x = get_locations(mesher, direction1)
  y = get_locations(mesher, direction2)

  dxMap = add!(mult!(FirstDerivativeOp(direction1, mesher), (-x * get_a(model))),
              mult!(SecondDerivativeOp(direction1, mesher), (0.5 * get_sigma(model) * get_sigma(model) * ones(mesher.layout.size))))

  dyMap = add!(mult!(FirstDerivativeOp(direction2, mesher), (-y * get_b(model))),
              mult!(SecondDerivativeOp(direction2, mesher), (0.5 * get_eta(model) * get_eta(model) * ones(mesher.layout.size))))

  corrMap = mult!(SecondOrderMixedDerivativeOp(direction1, direction2, mesher), fill(get_rho(model) * get_sigma(model) * get_eta(model), mesher.layout.size))

  mapX = TripleBandLinearOp(direction1, mesher)
  mapY = TripleBandLinearOp(direction2, mesher)

  return FdmG2Op(direction1, direction2, x, y, dxMap, dyMap, corrMap, mapX, mapY, model)
end

type FdmHullWhiteOp{I <: Integer} <: FdmLinearOpComposite
  direction::I
  x::Vector{Float64}
  dzMap::TripleBandLinearOp
  mapT::TripleBandLinearOp
  model::HullWhite
end

function FdmHullWhiteOp(mesher::FdmMesher, model::HullWhite, direction::Int)
  x = get_locations(mesher, direction)

  dzMap = add!(mult!(FirstDerivativeOp(direction, mesher), (-x * get_a(model))),
              mult!(SecondDerivativeOp(direction, mesher), (0.5 * get_sigma(model) * get_sigma(model)) * ones(mesher.layout.size)))

  mapT = TripleBandLinearOp(direction, mesher)

  return FdmHullWhiteOp(direction, x, dzMap, mapT, model)
end

# OP methods #
function set_time!(op::FdmG2Op, t1::Float64, t2::Float64)
  dynamics = get_dynamics(op.model)

  phi = 0.5 * (short_rate(dynamics, t1, 0.0, 0.0) + short_rate(dynamics, t2, 0.0, 0.0))

  hr = -0.5 * (op.x + op.y + phi)
  axpyb!(op.mapX, Vector{Float64}(0), op.dxMap, op.dxMap, hr)
  axpyb!(op.mapY, Vector{Float64}(0), op.dyMap, op.dyMap, hr)

  return op
end

function set_time!(op::FdmHullWhiteOp, t1::Float64, t2::Float64)
  dynamics = get_dynamics(op.model)

  phi = 0.5 * (short_rate(dynamics, t1, 0.0) + short_rate(dynamics, t2, 0.0))

  axpyb!(op.mapT, Vector{Float64}(0), op.dzMap, op.dzMap, -(op.x + phi))

  return op
end

function apply_direction(op::FdmG2Op, direction::Int, r::Vector{Float64})
  if direction == op.direction1
    return apply(op.mapX, r)
  elseif direction == op.direction2
    return apply(op.mapY, r)
  else
    return zeros(length(r))
  end
end

function apply_direction(op::FdmHullWhiteOp, direction::Int, r::Vector{Float64})
  if direction == op.direction
    return apply(op.mapT, r)
  else
    return zeros(length(r))
  end
end

function solve_splitting(op::FdmG2Op, direction::Int, r::Vector{Float64}, a::Float64)
  if direction == op.direction1
    return solve_splitting(op.mapX, r, a, 1.0)
  elseif direction == op.direction2
    return solve_splitting(op.mapY, r, a, 1.0)
  else
    return zeros(length(r))
  end
end

function solve_splitting(op::FdmHullWhiteOp, direction::Int, r::Vector{Float64}, a::Float64)
  if direction == op.direction
    return solve_splitting(op.mapT, r, a, 1.0)
  else
    return zeros(length(r))
  end
end

apply_mixed(op::FdmG2Op, r::Vector{Float64}) = apply(op.corrMap, r)
apply_mixed(op::FdmHullWhiteOp, r::Vector{Float64}) = zeros(length(r))

apply(op::FdmG2Op, r::Vector{Float64}) = apply(op.mapX, r) + apply(op.mapY, r) + apply_mixed(op, r)
apply(op::FdmHullWhiteOp, r::Vector{Float64}) = apply(op.mapT, r)

get_size(op::FdmG2Op) = 2
get_size(op::FdmHullWhiteOp) = 1

type FdmSolverDesc{F <: FdmMesher, C <: FdmInnerValueCalculator, I <: Integer}
  mesher::F
  bcSet::FdmBoundaryConditionSet
  condition::FdmStepConditionComposite
  calculator::C
  maturity::Float64
  timeSteps::I
  dampingSteps::I
end

## Schemes ##
type HundsdorferScheme <: FdScheme
  theta::Float64
  mu::Float64
  map::FdmLinearOpComposite
  bcSet::FdmBoundaryConditionSet
  dt::Float64
end

HundsdorferScheme(theta::Float64, mu::Float64, map::FdmLinearOpComposite, bcSet::FdmBoundaryConditionSet) = HundsdorferScheme(theta, mu, map, bcSet, 0.0)

type DouglasScheme <: FdScheme
  theta::Float64
  map::FdmLinearOpComposite
  bcSet::FdmBoundaryConditionSet
  dt::Float64
end

DouglasScheme(theta::Float64, map::FdmLinearOpComposite, bcSet::FdmBoundaryConditionSet) = DouglasScheme(theta, map, bcSet, 0.0)

set_step!(evolver::FdScheme, dt::Float64) = evolver.dt = dt

function step!(evolver::HundsdorferScheme, a::Vector{Float64}, t::Float64)
  t - evolver.dt > -1e-8 || error("a step towards negative time given")

  set_time!(evolver.map, max(0.0, t - evolver.dt), t)
  set_time!(evolver.bcSet, max(0.0, t - evolver.dt))

  apply_before_applying!(evolver.bcSet, evolver.map)
  y = a + evolver.dt * apply(evolver.map, a)
  apply_after_applying!(evolver.bcSet, y)

  y0 = copy(y)

  for i = 1:get_size(evolver.map)
    rhs = y - evolver.theta * evolver.dt * apply_direction(evolver.map, i, a)
    y = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_before_applying!(evolver.bcSet, evolver.map)
  yt = y0 + evolver.mu * evolver.dt * apply(evolver.map, y - a)
  apply_after_applying!(evolver.bcSet, yt)

  for i = 1:get_size(evolver.map)
    rhs = yt - evolver.theta * evolver.dt * apply_direction(evolver.map, i, y)
    yt = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_after_solving!(evolver.bcSet, yt)

  a[:] = yt

  return a, evolver
end

function step!(evolver::DouglasScheme, a::Vector{Float64}, t::Float64)
  t - evolver.dt > -1e-8 || error("a step towards negative time given")

  set_time!(evolver.map, max(0.0, t - evolver.dt), t)
  set_time!(evolver.bcSet, max(0.0, t - evolver.dt))

  apply_before_applying!(evolver.bcSet, evolver.map)
  y = a + evolver.dt * apply(evolver.map, a)
  apply_after_applying!(evolver.bcSet, y)

  for i = 1:get_size(evolver.map)
    rhs = y - evolver.theta * evolver.dt * apply_direction(evolver.map, i, a)
    y = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_after_solving!(evolver.bcSet, y)

  a[:] = y

  return a, evolver
end

## Main Finite Difference Model ##
type FiniteDifferenceModel{T}
  evolver::T
  stoppingTimes::Vector{Float64}

  FiniteDifferenceModel{T}(evolver::T, stoppingTimes::Vector{Float64}) = new(evolver, sort(unique(stoppingTimes)))
end

function rollback_impl!(model::FiniteDifferenceModel, a::Vector{Float64}, from::Float64, to::Float64, steps::Int, condition::StepCondition)
  dt = (from - to) / steps
  t = from
  set_step!(model.evolver, dt)

  if length(model.stoppingTimes) > 0 && model.stoppingTimes[end] == from
    apply_to!(condition, a, from)
  end

  for i = 1:steps
    _now = t
    _next = t - dt
    hit = false
    for j = length(model.stoppingTimes):-1:1
      if _next <= model.stoppingTimes[j] && model.stoppingTimes[j] < _now
        # a stopping time was hit
        hit = true

        # perform a small step to stoppingTimes[j]
        set_step!(model.evolver, _now - model.stoppingTimes[j])
        step!(model.evolver, a, _now)

        apply_to!(condition, a, model.stoppingTimes[j])

        # and continue the cycle
        _now = model.stoppingTimes[j]
      end
    end

    # if we did hit
    if hit
      # we might have to make a small step to complete the big one
      if _now > _next
        set_step!(model.evolver, _now - _next)
        step!(model.evolver, a, _now)
        apply_to!(condition, a, _next)
      end

      # and in any case, we have to reset the evolver to the default step
      set_step!(model.evolver, dt)
    else
      # if we didn't, the evolver is already set to the default step, which is ok
      step!(model.evolver, a, _now)
      apply_to!(condition, a, _next)
    end
    t -= dt
  end

  return model, a
end

rollback!(model::FiniteDifferenceModel, a::Vector{Float64}, from::Float64, to::Float64, steps::Int, condition::StepCondition) = rollback_impl!(model, a, from, to, steps, condition)

## Solvers ##
type FdmBackwardSolver
  map::FdmLinearOpComposite
  bcSet::FdmBoundaryConditionSet
  condition::FdmStepConditionComposite
  schemeDesc::FdmSchemeDesc
end

function FdmBackwardSolver(map::FdmLinearOpComposite, bcSet::FdmBoundaryConditionSet, schemeDesc::FdmSchemeDesc)
  condition = FdmStepConditionComposite(Vector{Float64}(), Vector{StepCondition}())

  return FdmBackwardSolver(map, bcSet, condition, schemeDesc)
end

function rollback!(bsolv::FdmBackwardSolver, schemeType::Hundsdorfer, rhs::Vector{Float64}, from::Float64, to::Float64, steps::Int, dampingSteps::Int)
  deltaT = from - to
  allSteps = steps + dampingSteps
  dampingTo = from - (deltaT * dampingSteps) / allSteps

  # if dampingSteps > 0
  #   implicitEvolver = ImplicitEulerScheme(bsolv.map, bsolv.bcSet)
  # end

  dampingSteps == 0 || error("damping steps shoudl be 0")

  hsEvolver = HundsdorferScheme(bsolv.schemeDesc.theta, bsolv.schemeDesc.mu, bsolv.map, bsolv.bcSet)
  hsModel = FiniteDifferenceModel{HundsdorferScheme}(hsEvolver, bsolv.condition.stoppingTimes)
  rollback!(hsModel, rhs, dampingTo, to, steps, bsolv.condition)

  return bsolv, rhs
end

function rollback!(bsolv::FdmBackwardSolver, schemeType::Douglas, rhs::Vector{Float64}, from::Float64, to::Float64, steps::Int, dampingSteps::Int)
  deltaT = from - to
  allSteps = steps + dampingSteps
  dampingTo = from - (deltaT * dampingSteps) / allSteps

  dampingSteps == 0 || error("damping steps shoudl be 0")

  dsEvolver = DouglasScheme(bsolv.schemeDesc.theta, bsolv.map, bsolv.bcSet)
  dsModel = FiniteDifferenceModel{DouglasScheme}(dsEvolver, bsolv.condition.stoppingTimes)
  rollback!(dsModel, rhs, dampingTo, to, steps, bsolv.condition)

  return bsolv, rhs
end

type Fdm1DimSolver{FD <: FdmLinearOpComposite} <: LazyObject
  lazyMixin::LazyMixin
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  op::FD
  thetaCondition::FdmSnapshotCondition
  conditions::FdmStepConditionComposite
  initialValues::Vector{Float64}
  resultValues::Vector{Float64}
  x::Vector{Float64}
  interpolation::NaturalCubicSpline

  function Fdm1DimSolver{FD}(solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc, op::FD)
    thetaCondition = FdmSnapshotCondition(0.99 * min(1.0 / 365.0, length(solverDesc.condition.stoppingTimes) == 0 ? solverDesc.maturity : solverDesc.condition.stoppingTimes[1]))
    conditions = join_conditions_FdmStepConditionComposite(thetaCondition, solverDesc.condition)

    layout = solverDesc.mesher.layout

    x = zeros(layout.size)
    initialValues = zeros(layout.size)
    resultValues = zeros(layout.size)

    coords = ones(Int, length(layout.dim))

    for i = 1:layout.size
      initialValues[i] = avg_inner_value(solverDesc.calculator, coords, i, solverDesc.maturity)
      x[i] = get_location(solverDesc.mesher, coords, 1)

      iter_coords!(coords, layout.dim)
    end

    new(LazyMixin(), solverDesc, schemeDesc, op, thetaCondition, conditions, initialValues, resultValues, x)
  end
end

type Fdm2DimSolver{FD <: FdmLinearOpComposite} <: LazyObject
  lazyMixin::LazyMixin
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  op::FD
  thetaCondition::FdmSnapshotCondition
  conditions::FdmStepConditionComposite
  initialValues::Vector{Float64}
  resultValues::Matrix{Float64}
  x::Vector{Float64}
  y::Vector{Float64}
  interpolation::BicubicSpline

  function Fdm2DimSolver{FD}(solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc, op::FD)
    thetaCondition = FdmSnapshotCondition(0.99 * min(1.0 / 365.0, length(solverDesc.condition.stoppingTimes) == 0 ? solverDesc.maturity : solverDesc.condition.stoppingTimes[1]))
    conditions = join_conditions_FdmStepConditionComposite(thetaCondition, solverDesc.condition)

    layout = solverDesc.mesher.layout

    initialValues = zeros(layout.size)
    resultValues = zeros(layout.dim[1], layout.dim[2])

    x = zeros(layout.dim[1])
    y = zeros(layout.dim[2])
    x_count = 1
    y_count = 1

    coords = ones(Int, length(layout.dim))

    for i = 1:layout.size
      initialValues[i] = avg_inner_value(solverDesc.calculator, coords, i, solverDesc.maturity)

      if coords[2] == 1
        x[x_count] = get_location(solverDesc.mesher, coords, 1)
        x_count += 1
      end

      if coords[1] == 1
        y[y_count] =  get_location(solverDesc.mesher, coords, 2)
        y_count += 1
      end

      iter_coords!(coords, layout.dim)
    end

    return new{FD}(LazyMixin(), solverDesc, schemeDesc, op, thetaCondition, conditions, initialValues, resultValues, x, y)
  end
end

get_interpolation(solv::Fdm1DimSolver, x::Float64) = solv.interpolation(x)
get_interpolation(solv::Fdm2DimSolver, x::Float64, y::Float64) = solv.interpolation.spline(x, y)

function interpolate_at(solv::Fdm1DimSolver, x::Float64)
  calculate!(solv)
  return get_interpolation(solv, x)
end

function interpolate_at(solv::Fdm2DimSolver, x::Float64, y::Float64)
  calculate!(solv)
  return get_interpolation(solv, x, y)
end

function perform_calculations!(solv::Fdm1DimSolver)
  rhs = copy(solv.initialValues)

  bsolver = FdmBackwardSolver(solv.op, solv.solverDesc.bcSet, solv.conditions, solv.schemeDesc)
  rollback!(bsolver, solv.schemeDesc.schemeType, rhs, solv.solverDesc.maturity, 0.0, solv.solverDesc.timeSteps, solv.solverDesc.dampingSteps)
  solv.resultValues = copy(rhs)
  solv.interpolation = NaturalCubicSpline(solv.x, solv.resultValues)

  return solv
end

function perform_calculations!(solv::Fdm2DimSolver)
  rhs = copy(solv.initialValues)

  bsolver = FdmBackwardSolver(solv.op, solv.solverDesc.bcSet, solv.conditions, solv.schemeDesc)
  rollback!(bsolver, solv.schemeDesc.schemeType, rhs, solv.solverDesc.maturity, 0.0, solv.solverDesc.timeSteps, solv.solverDesc.dampingSteps)
  solv.resultValues = reshape(rhs, length(solv.x), length(solv.y))
  solv.interpolation = BicubicSpline(solv.x, solv.y, solv.resultValues)

  return solv
end

type FdmG2Solver <: LazyObject
  lazyMixin::LazyMixin
  model::G2
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  solver::Fdm2DimSolver

  function FdmG2Solver(model::G2, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc)
    op = FdmG2Op(solverDesc.mesher, model, 1, 2)
    solver = Fdm2DimSolver{FdmG2Op}(solverDesc, schemeDesc, op)
    new(LazyMixin(), model, solverDesc, schemeDesc, solver)
  end
end

type FdmHullWhiteSolver <: LazyObject
  lazyMixin::LazyMixin
  model::HullWhite
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  solver::Fdm1DimSolver

  function FdmHullWhiteSolver(model::HullWhite, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc)
    op = FdmHullWhiteOp(solverDesc.mesher, model, 1)
    solver = Fdm1DimSolver{FdmHullWhiteOp}(solverDesc, schemeDesc, op)
    new(LazyMixin(), model, solverDesc, schemeDesc, solver)
  end
end

function value_at(solv::FdmG2Solver, x::Float64, y::Float64)
  return interpolate_at(solv.solver, x, y)
end

function value_at(solv::FdmHullWhiteSolver, x::Float64)
  return interpolate_at(solv.solver, x)
end
