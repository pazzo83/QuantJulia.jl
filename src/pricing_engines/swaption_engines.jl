## Swaption Pricing Engines ##
using Distributions
using QuantJulia.Time, QuantJulia.Math

type BlackSwaptionEngine{Y <: YieldTermStructure, S <: SwaptionVolatilityStructure, DC <: DayCount} <: PricingEngine{Y}
  yts::Y
  vol::Quote
  volStructure::S
  dc::DC
  displacement::Float64
end

BlackSwaptionEngine{Y <: YieldTermStructure, DC <: DayCount}(yts::Y, vol::Quote, dc::DC, displacement::Float64 = 0.0) =
                    BlackSwaptionEngine(yts, vol, ConstantSwaptionVolatility(0, QuantJulia.Time.NullCalendar(), QuantJulia.Time.Following(), vol, dc), dc, displacement)

type G2SwaptionEngine{Y <: YieldTermStructure, I <: Integer} <: PricingEngine{Y}
  model::G2{Y}
  range::Float64
  intervals::I
end

type JamshidianSwaptionEngine{S <: ShortRateModel, Y <: YieldTermStructure} <: PricingEngine{Y}
  model::S
  ts::Y

  function call{S}(::Type{JamshidianSwaptionEngine}, model::S)
    new{S, YieldTermStructure}(model)
  end

  function call{S, Y}(::Type{JamshidianSwaptionEngine}, model::S, yts::Y)
    new{S, Y}(model, yts)
  end
end

type TreeSwaptionEngine{S <: ShortRateModel, I <: Integer, T <: ShortRateTree, Y <: YieldTermStructure} <: LatticeShortRateModelEngine{S, Y, T}
  model::S
  timeSteps::I
  tg::TimeGrid
  lattice::T
  ts::Y

  call{S, I, T}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, tg::TimeGrid, l::T) = new{S, I, T, YieldTermStructure}(m, tsteps, tg, l)

  call{S, I, T, Y}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, tg::TimeGrid, l::T, ts::Y) = new{S, I, T, Y}(m, tsteps, tg, l, ts)
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, tg::TimeGrid)
  lattice = tree(model, tg)
  ts = TreeSwaptionEngine(model, 0, tg, lattice)

  add_observer!(model, ts)

  return ts
end

function update!(eng::LatticeShortRateModelEngine)
  if length(eng.tg.times) > 0
    eng.lattice = tree(eng.model, eng.tg)
  end

  return eng
end

## General Pricing Functions ##
function black_formula{T <: OptionType}(optionType::T, strike::Float64, forward::Float64, stdDev::Float64, discount::Float64 = 1.0, displacement::Float64 = 0.0)
  # TODO check requirements (see cpp)
  opt_type = QuantJulia.value(optionType)
  if stdDev == 0.0
    return max((forward - strike) * opt_type, 0.0 * discount)
  end

  forward += displacement
  strike += displacement

  if strike == 0.0
    return isa(optionType, Call) ? forward * discount : 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev
  d2 = d1 - stdDev
  norm = Normal() # using distributions.jl
  nd1 = cdf(norm, opt_type * d1)
  nd2 = cdf(norm, opt_type * d2)
  result = discount * opt_type * (forward * nd1 - strike * nd2)

  return result
end

function black_formula_standard_dev_derivative(strike::Float64, forward::Float64, stdDev::Float64, discount::Float64, displacement::Float64)
  forward += displacement
  strike += displacement

  if stdDev == 0.0 || strike == 0.0
    return 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev

  return discount * forward * distribution_derivative(Normal(), d1)
end

get_annuity(delivery::SettlementPhysical, swap::VanillaSwap) = abs(fixed_leg_BPS(swap)) / basisPoint

function _calculate!(pe::BlackSwaptionEngine, swaption::Swaption)
  exerciseDate = swaption.exercise.dates[1]
  swap = swaption.swap

  strike = swap.fixedRate

  # override swap's pricing engine temporarily, bypassing normal calc flow, since swap.iborIndex might be using a diff curve
  _calculate!(DiscountingSwapEngine(pe.yts), swap)
  swap.lazyMixin.calculated = true
  atmForward = fair_rate(swap)

  if swap.spread != 0.0
    correction = swap.spread * abs(floating_leg_BPS(swap) / fixed_leg_BPS(swap))
    strike -= correction
    atmForward -= correction
    swaption.results.additionalResults["spreadCorrection"] = correction
  else
    swaption.results.additionalResults["spreadCorrection"] = 0.0
  end

  swaption.results.additionalResults["strike"] = strike
  swaption.results.additionalResults["atmForward"] = atmForward

  annuity = get_annuity(swaption.delivery, swap)
  swaption.results.additionalResults["annuity"] = annuity

  # the swap length calculation might be improved using the value date of the exercise date
  swapLength = swap_length(pe.volStructure, exerciseDate, date(get_latest_coupon(swap.legs[1])))
  swaption.results.additionalResults["swapLength"] = swapLength

  variance = black_varience(pe.volStructure, exerciseDate, swapLength, strike)

  stdDev = sqrt(variance)
  swaption.results.additionalResults["stdDev"] = stdDev
  w = isa(swap.swapT, Payer) ? Call() : Put()

  swaption.results.value = black_formula(w, strike, atmForward, stdDev, annuity, pe.displacement)

  exerciseTime = time_from_reference(pe.volStructure, exerciseDate)

  swaption.results.additionalResults["vega"] = sqrt(exerciseTime) * black_formula_standard_dev_derivative(strike, atmForward, stdDev, annuity, pe.displacement)

  # resetting swap
  reset!(swap.results)
  swap.lazyMixin.calculated = false

  return swaption
end

function _calculate!(pe::G2SwaptionEngine, swaption::Swaption)
  swap = swaption.swap

  # overriding pricing engine
  _calculate!(DiscountingSwapEngine(pe.model.ts), swap)
  swap.lazyMixin.calculated = true

  correction = swap.spread * abs(floating_leg_BPS(swap) / fixed_leg_BPS(swap))
  fixedRate = swap.fixedRate - correction
  swaption.results.value = gen_swaption(pe.model, swaption, fixedRate, pe.range, pe.intervals)

  return swaption
end

function _calculate!(pe::JamshidianSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model
  ref_date = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  amounts = copy(swaption.swap.args.fixedCoupons)
  amounts[end] += swaption.swap.nominal

  maturity = year_fraction(dc, ref_date, swaption.exercise.dates[1])

  fixedPayTimes = zeros(length(swaption.swap.args.fixedPayDates))
  valueTime = year_fraction(dc, ref_date, swaption.swap.args.fixedResetDates[1])
  for i = 1:length(fixedPayTimes)
    fixedPayTimes[i] = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])
  end

  finder = RStarFinder(tsmodel, swaption.swap.nominal, maturity, valueTime, fixedPayTimes, amounts)

  minStrike = -10.0
  maxStrike = 10.0
  slv = BrentSolver(10000, true, true, minStrike, maxStrike)

  rStar = solve(slv, operator(finder), 1e-8, 0.05, minStrike, maxStrike)

  w = isa(swaption.swap.swapT, Payer) ? Put() : Call()

  _size = length(swaption.swap.args.fixedCoupons)

  val = 0.0
  _B = discount_bond(tsmodel, maturity, valueTime, rStar)

  @simd for i = 1:_size
    @inbounds fixedPayTime = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])

    strike = discount_bond(tsmodel, maturity, fixedPayTime, rStar) / _B

    dboValue = discount_bond_option(tsmodel, w, strike, maturity, valueTime, fixedPayTime)

    @inbounds val += amounts[i] * dboValue
  end

  swaption.results.value = val

  return swaption
end

function greater_than_or_equal_to{T}(x::T, y::T)
  return x >= y
end

function _calculate!(pe::TreeSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  dSwaption = DiscretizedSwaption(swaption, refDate, dc)

  lattice = pe.lattice
  stoppingTimes = zeros(length(swaption.exercise.dates))
  @simd for i = 1:length(stoppingTimes)
    @inbounds stoppingTimes[i] = year_fraction(dc, refDate, swaption.exercise.dates[i])
  end

  initialize!(dSwaption, lattice.treeLattice, stoppingTimes[end])

  nextExerciseIdx = findnext(greater_than_or_equal_to, stoppingTimes, 1, 0.0)

  nextExercise = nextExerciseIdx != 0 ? stoppingTimes[nextExerciseIdx] : stoppingTimes[end]

  rollback!(dSwaption, nextExercise)
  swaption.results.value = present_value(dSwaption)
end
