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

get_dminus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dminus[coords[direction]]
get_dplus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dplus[coords[direction]]

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


type FdmSolverDesc{F <: FdmMesher, C <: FdmInnerValueCalculator, I <: Integer}
  mesher::F
  bcSet::FdmBoundaryConditionSet
  condition::FdmStepConditionComposite
  calculator::C
  maturity::Float64
  timeSteps::I
  dampingSteps::I
end

# Not done
type Fdm2DimSolver{FD <: FdmLinearOpComposite}
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  op::FD
end

type FdmG2Solver <: LazyObject
  lazyMixin::LazyMixin
  model::G2
  solverDesc::FdmSolverDesc
  schemeDesc::FdmSchemeDesc
  solver::Fdm2DimSolver

  function FdmG2Solver(model::G2, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc)
    op = FdmG2Op(solverDesc.mesher, model, 1, 2)
    solver = Fdm2DimSolver(solverDesc, schemeDesc, op)
    new(LazyMixin(), model, solverDesc, schemeDesc, solver)
  end
end
