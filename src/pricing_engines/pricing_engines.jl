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
  model::ShortRateModel
  ts::Y

  function call{S}(::Type{JamshidianSwaptionEngine}, model::S)
    new{S, YieldTermStructure}(model)
  end

  function call{S, Y}(::Type{JamshidianSwaptionEngine}, model::S, yts::Y)
    new{S, Y}(model, yts)
  end
end

type DiscretizedSwap <: DiscretizedAsset
  fixedResetTimes::Vector{Float64}
  fixedPayTimes::Vector{Float64}
  floatingResetTimes::Vector{Float64}
  floatingPayTimes::Vector{Float64}

  function DiscretizedSwap{DC <: DayCount}(referenceDate::Date, dc::DC, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date}, floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date})
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

    new(fixedResetTimes, fixedPayTimes, floatingResetTimes, floatingPayTimes)
  end
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
end

function DiscretizedSwaption{DC <: DayCount}(swaption::Swaption, referenceDate::Date, dc::DC)
  dates = swaption.exercise.dates
  fixed_coups = swaption.swap.legs[1].coupons
  floating_coups = swaption.swap.legs[2].coupons
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
  for i = 1:n
    exerciseDate = dates[i]

    for j = 1:length(fixed_coups)
      if within_next_week(exerciseDate, fixedPayDates[j]) && fixedResetDates[j] < referenceDate
        fixedPayDates[j] = exerciseDate
      end

      if within_previous_week(exerciseDate, fixedResetDates[j])
        fixedResetDates[j] = exerciseDate
      end
    end

    for j = 1:length(floating_coups)
      if within_previous_week(exerciseDate, floatingResetDates[i])
        floatingResetDates[j] = exerciseDate
      end
    end
  end

  lastFixedPayment = year_fraction(dc, referenceDate, fixedPayDates[end])
  lastFloatingPayment = year_fraction(dc, referenceDate, floatingPayDates[end])

  lastPayment = max(lastFixedPayment, lastFloatingPayment)
  underlying = DiscretizedSwap(referenceDate, dc, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates)
  exercise = swaption.exercise

  DiscretizedSwaption(underlying, exercise, exerciseTimes, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, lastPayment)
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

  for i = 1:_size
    fixedPayTime = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])

    strike = discount_bond(tsmodel, maturity, fixedPayTime, rStar) / _B

    dboValue = discount_bond_option(tsmodel, w, strike, maturity, valueTime, fixedPayTime)

    val += amounts[i] * dboValue
  end

  swaption.results.value = val

  return swaption
end
