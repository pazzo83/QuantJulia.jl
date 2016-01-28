# Pricing Engines module
# module PricingEngines
using Distributions
using QuantJulia.Time, QuantJulia.Math

const basisPoint = 0.0001

type DiscountingBondEngine{Y <: YieldTermStructure} <: PricingEngine{Y}
  yts::Y

  function call(::Type{DiscountingBondEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingBondEngine}, yts::Y)
    new{Y}(yts)
  end
end

type DiscountingSwapEngine{Y <: YieldTermStructure} <: PricingEngine{Y}
  yts::Y

  function call(::Type{DiscountingSwapEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingSwapEngine}, yts::Y)
    new{Y}(yts)
  end
end

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

type DiscretizedAssetCommon{L <: Lattice}
  time::Float64
  values::Vector{Float64}
  latestPreAdjustment::Float64
  latestPostAdjustment::Float64
  method::L

  call(::Type{DiscretizedAssetCommon}, t::Float64, v::Vector{Float64}, lpa::Float64, lpoa::Float64) =
      new{Lattice}(t, v, lpa, lpoa)
end

DiscretizedAssetCommon() = DiscretizedAssetCommon(0.0, zeros(0), eps(), eps())

set_time!(a::DiscretizedAsset, t::Float64) = a.common.time = t
set_method!(a::DiscretizedAsset, method::Lattice) = a.common.method = method

type DiscretizedSwap{ST <: SwapType} <: DiscretizedAsset
  nominal::Float64
  swapT::ST
  fixedResetTimes::Vector{Float64}
  fixedPayTimes::Vector{Float64}
  floatingResetTimes::Vector{Float64}
  floatingPayTimes::Vector{Float64}
  args::VanillaSwapArgs
  common::DiscretizedAssetCommon
end

function DiscretizedSwap{DC <: DayCount, ST <: SwapType}(nominal::Float64, swapT::ST, referenceDate::Date, dc::DC, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date}, floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date}, args::VanillaSwapArgs)
  fixed_n = length(fixedPayDates)
  float_n = length(floatingPayDates)

  fixedResetTimes = zeros(fixed_n)
  fixedPayTimes = zeros(fixed_n)
  floatingResetTimes = zeros(float_n)
  floatingPayTimes = zeros(float_n)

  for i = 1:fixed_n
    fixedResetTimes[i] = year_fraction(dc, referenceDate, fixedResetDates[i])
    fixedPayTimes[i] = year_fraction(dc, referenceDate, fixedPayDates[i])
  end

  for i = 1:float_n
    floatingResetTimes[i] = year_fraction(dc, referenceDate, floatingResetDates[i])
    floatingPayTimes[i] = year_fraction(dc, referenceDate, floatingPayDates[i])
  end

  DiscretizedSwap(nominal, swapT, fixedResetTimes, fixedPayTimes, floatingResetTimes, floatingPayTimes, args, DiscretizedAssetCommon())
end

type DiscretizedSwaption{E <: Exercise} <: DiscretizedOption
  underlying::DiscretizedSwap
  exercise::E
  exerciseTimes::Vector{Float64}
  fixedPayDates::Vector{Date}
  fixedResetDates::Vector{Date}
  floatingPayDates::Vector{Date}
  floatingResetDates::Vector{Date}
  lastPayment::Float64
  common::DiscretizedAssetCommon
end

function DiscretizedSwaption{DC <: DayCount}(swaption::Swaption, referenceDate::Date, dc::DC)
  dates = swaption.exercise.dates
  fixed_coups = swaption.swap.legs[1].coupons
  floating_coups = swaption.swap.legs[2].coupons
  nominal = swaption.swap.nominal
  swapT = swaption.swap.swapT
  n = length(dates)

  exerciseTimes = zeros(n)
  # fixedPayDates = get_pay_dates(fixed_coups)
  # fixedResetDates = get_reset_dates(fixed_coups)
  # floatingPayDates = get_pay_dates(floating_coups)
  # floatingResetDates = get_reset_dates(floating_coups)
  fixedPayDates = swaption.swap.args.fixedPayDates
  fixedResetDates = swaption.swap.args.fixedResetDates
  floatingPayDates = swaption.swap.args.floatingPayDates
  floatingResetDates = swaption.swap.args.floatingResetDates

  for i = 1:n
    exerciseTimes[i] = year_fraction(dc, referenceDate, dates[i])
  end

  # Date adjustments can get time vectors out of sync
  # Here we try and collapse similar dates which could cause a mispricing
  @simd for i = 1:n
    @inbounds exerciseDate = dates[i]

    for j = 1:length(fixed_coups)
      @inbounds if within_next_week(exerciseDate, fixedPayDates[j]) && fixedResetDates[j] < referenceDate
        @inbounds fixedPayDates[j] = exerciseDate
      end

      @inbounds if within_previous_week(exerciseDate, fixedResetDates[j])
        @inbounds fixedResetDates[j] = exerciseDate
      end
    end

    for j = 1:length(floating_coups)
      @inbounds if within_previous_week(exerciseDate, floatingResetDates[i])
        @inbounds floatingResetDates[j] = exerciseDate
      end
    end
  end

  lastFixedPayment = year_fraction(dc, referenceDate, fixedPayDates[end])
  lastFloatingPayment = year_fraction(dc, referenceDate, floatingPayDates[end])

  lastPayment = max(lastFixedPayment, lastFloatingPayment)
  underlying = DiscretizedSwap(nominal, swapT, referenceDate, dc, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, swaption.swap.args)
  exercise = swaption.exercise

  DiscretizedSwaption(underlying, exercise, exerciseTimes, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, lastPayment, DiscretizedAssetCommon())
end

type DiscretizedDiscountBond <: DiscretizedAsset
  common::DiscretizedAssetCommon
end

DiscretizedDiscountBond() = DiscretizedDiscountBond(DiscretizedAssetCommon())

function adjust_values!(dAsset::DiscretizedAsset)
  pre_adjust_values!(dAsset)
  post_adjust_values!(dAsset)

  return dAsset
end

function pre_adjust_values!(dAsset::DiscretizedAsset)
  if ~QuantJulia.Math.close_enough(dAsset.common.time, dAsset.common.latestPreAdjustment)
    pre_adjust_values_impl!(dAsset)
    dAsset.common.latestPreAdjustment = dAsset.common.time
  end

  return dAsset
end

function post_adjust_values!(dAsset::DiscretizedAsset)
  if ~QuantJulia.Math.close_enough(dAsset.common.time, dAsset.common.latestPostAdjustment)
    post_adjust_values_impl!(dAsset)
    dAsset.common.latestPostAdjustment = dAsset.common.time
  end

  return dAsset
end

function partial_rollback!(dAsset::DiscretizedAsset, t::Float64)
  partial_rollback!(dAsset.common.method, dAsset, t)
  return dAsset
end

function is_on_time(dAsset::DiscretizedAsset, t::Float64)
  grid = dAsset.common.method.tg
  return QuantJulia.Math.close_enough(grid.times[findfirst(grid.times .>= t)], dAsset.common.time)
end

function rollback!(dAsset::DiscretizedAsset, t::Float64)
  rollback!(dAsset.common.method, dAsset, t)

  return dAsset
end

pre_adjust_values_impl!(dAsset::DiscretizedAsset) = dAsset # do nothing
post_adjust_values_impl!(dAsset::DiscretizedAsset) = dAsset # do nothing

function general_reset!(dOption::DiscretizedOption, sz::Int)
  # check methods in main and underlying
  dOption.common.values = zeros(sz)
  adjust_values!(dOption)
  return dOption
end

present_value(dAsset::DiscretizedAsset) = present_value(dAsset.common.method, dAsset)


function post_adjust_values_impl!(dOption::DiscretizedOption)
  partial_rollback!(dOption.underlying, dOption.common.time)
  pre_adjust_values!(dOption.underlying)

  if isa(dOption.exercise, AmericanExercise)
    if dOption.common.time >= dOption.exerciseTimes[1] && dOption.common.time <= dOption.exerciseTimes[1]
      apply_exercise_condition!(dOption)
    end
  else
    @simd for i = 1:length(dOption.exerciseTimes)
      @inbounds t = dOption.exerciseTimes[i]
      if t >= 0.0 && is_on_time(dOption, t)
        apply_exercise_condition!(dOption)
      end
    end
  end

  post_adjust_values!(dOption.underlying)

  return dOption
end

function apply_exercise_condition!(dOption::DiscretizedOption)
  @simd for i = 1:length(dOption.common.values)
    @inbounds dOption.common.values[i] = max(dOption.underlying.common.values[i], dOption.common.values[i])
  end

  return dOption
end

function mandatory_times(discretizedSwap::DiscretizedSwap)
  # get times
  times = vcat(discretizedSwap.fixedResetTimes[discretizedSwap.fixedResetTimes .>= 0.0], discretizedSwap.fixedPayTimes[discretizedSwap.fixedPayTimes .>= 0.0],
               discretizedSwap.floatingResetTimes[discretizedSwap.floatingResetTimes .>= 0.0], discretizedSwap.floatingPayTimes[discretizedSwap.floatingPayTimes .>= 0.0])

  return times
end

function mandatory_times(discretizedSwaption::DiscretizedSwaption)
  times = mandatory_times(discretizedSwaption.underlying)
  times =  times[times .>= 0.0]
  times = vcat(times, discretizedSwaption.exerciseTimes)
  return times
end

function initialize!(asset::DiscretizedAsset, lattice::Lattice, t::Float64)
  set_method!(asset, lattice)
  initialize!(lattice, asset, t)
  return asset
end

function reset!(dSwaption::DiscretizedSwaption, sz::Int)
  initialize!(dSwaption.underlying, dSwaption.common.method, dSwaption.lastPayment)
  general_reset!(dSwaption, sz)
  return dSwaption
end

function reset!(dSwap::DiscretizedSwap, sz::Int)
  dSwap.common.values = zeros(sz)
  adjust_values!(dSwap)

  return dSwap
end

function reset!(dBond::DiscretizedDiscountBond, sz::Int)
  dBond.common.values = ones(sz)

  return dBond
end

function pre_adjust_values_impl!(dSwap::DiscretizedSwap)
  # Floating payments
  @simd for i = 1:length(dSwap.floatingResetTimes)
    @inbounds t = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond()
      initialize!(bond, dSwap.common.method, dSwap.floatingPayTimes[i])
      rollback!(bond, dSwap.common.time)

      nominal = dSwap.nominal
      T = dSwap.args.floatingAccrualTimes[i]
      spread = dSwap.args.floatingSpreads[i]

      accruedSpread = nominal * T * spread
      for j = 1:length(dSwap.common.values)
        @inbounds coup = nominal * (1.0 - bond.common.values[j]) + accruedSpread * bond.common.values[j]

        if isa(dSwap.swapT, Payer)
          @inbounds dSwap.common.values[j] += coup
        else
          @inbounds dSwap.common.values[j] -= coup
        end
      end
    end
  end

  # Fixed Payments
  @simd for i = 1:length(dSwap.fixedResetTimes)
    t = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond()
      initialize!(bond, dSwap.common.method, dSwap.fixedPayTimes[i])
      rollback!(bond, dSwap.common.time)

      @inbounds fixedCoup = dSwap.args.fixedCoupons[i]

      for j = 1:length(dSwap.common.values)
        @inbounds coup = fixedCoup * bond.common.values[j]
        if isa(dSwap.swapT, Payer)
          @inbounds dSwap.common.values[j] -= coup
        else
          @inbounds dSwap.common.values[j] += coup
        end
      end
    end
  end

  return dSwap
end

function post_adjust_values_impl!(dSwap::DiscretizedSwap)
  # fixed coupons whose reset time is in the past won't be managed in pre_adjust_values
  @simd for i = 1:length(dSwap.fixedPayTimes)
    @inbounds t = dSwap.fixedPayTimes[i]
    @inbounds _reset = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      fixedCoup = dSwap.args.fixedCoupons[i]
      if isa(dSwap.swapT, Payer)
        dSwap.common.values -= fixedCoup
      else
        dSwap.common.values += fixedCoup
      end
    end
  end

  # the same applies to floating payments whose rate is already fixed
  @simd for i = 1:length(dSwap.floatingPayTimes)
    @inbounds t = dSwap.floatingPayTimes[i]
    @inbounds _reset = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      @inbounds currentFloatingCoup = dSwap.args.floatingCoupons[i]

      if isa(dSwap.swapT, Payer)
        dSwap.common.values += currentFloatingCoup
      else
        dSwap.common.values -= currentFloatingCoup
      end
    end
  end

  return dSwap
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


function _calculate!{B <: Bond}(pe::DiscountingBondEngine, bond::B)
  yts = pe.yts
  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

function _calculate!{S <: Swap}(pe::DiscountingSwapEngine, swap::S)
  # stuff
  # println("NEW ONE=============================================================================")
  # if swap.rate.value > 0.0323
  #   error("DUIE")
  # end
  swap.results.value = 0.0
  yts = pe.yts

  ref_date = yts.referenceDate

  swap.results.npvDateDiscount = discount(yts, ref_date)

  # for (i, leg) in enumerate(swap.legs)
  for i = 1:length(swap.legs)
    leg = swap.legs[i]
    swap.results.legNPV[i], swap.results.legBPS[i] = npvbps(leg, yts, ref_date, ref_date)
    swap.results.legNPV[i] *= swap.payer[i]
    swap.results.legBPS[i] *= swap.payer[i]

    d1 = accrual_start_date(leg.coupons[1])
    if d1 >= ref_date
      swap.results.startDiscounts[i] = discount(yts, d1)
    end

    d2 = accrual_end_date(leg.coupons[end])
    if (d2 >= ref_date)
      swap.results.endDiscounts[i] = discount(yts, d2)
    end

    swap.results.value += swap.results.legNPV[i]
  end

  if swap.results.legBPS[1] != 0.0
    swap.results.fairRate = swap.fixedRate - swap.results.value / (swap.results.legBPS[1] / basisPoint)
  end

  return swap
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
