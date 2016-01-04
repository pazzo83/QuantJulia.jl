# cash_flows.jl
# module CF

using QuantJulia.Time

## TYPES ##
## Coupon pricers
type BlackIborCouponPricer <: IborCouponPricer
  discount::Float64
  spreadLegValue::Float64
  initialized::Bool
end

BlackIborCouponPricer() = BlackIborCouponPricer(0.0, 0.0, false)


## types of cash flows ##
type SimpleCashFlow <: CashFlow
  amount::Float64
  date::Date
end

amount(cf::SimpleCashFlow) = cf.amount
date(cf::SimpleCashFlow) = cf.date

type FixedRateCoupon{DC <: DayCount} <: Coupon
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
  accrualStartDate::Date
  accrualEndDate::Date
  refPeriodStart::Date
  refPeriodEnd::Date
  dc::DC
end

type IborCoupon{I <: Integer, DC <: DayCount, ICP <: IborCouponPricer} <: Coupon
  paymentDate::Date
  nominal::Float64
  accrualStartDate::Date
  accrualEndDate::Date
  fixingDate::Date
  fixingValueDate::Date
  fixingEndDate::Date
  fixingDays::I
  iborIndex::IborIndex
  gearing::Float64
  spread::Float64
  refPeriodStart::Date
  refPeriodEnd::Date
  dc::DC
  isInArrears::Bool
  spanningTime::Float64
  pricer::ICP
end

function IborCoupon{I <: Integer, DC <: DayCount, ICP <: IborCouponPricer}(paymentDate::Date, nominal::Float64, startDate::Date, endDate::Date, fixingDays::I, iborIndex::IborIndex,
                    gearing::Float64, spread::Float64, refPeriodStart::Date, refPeriodEnd::Date, dc::DC, isInArrears::Bool,
                    pricer::ICP)
  # TODO check if right fixing days
  _fixing_date = fixing_date(iborIndex, startDate)
  fixing_cal = iborIndex.fixingCalendar
  idx_fixing_days = iborIndex.fixingDays
  fixing_val_date = advance(Base.Dates.Day(idx_fixing_days), fixing_cal, _fixing_date, iborIndex.convention)

  if isInArrears
    fixing_end_date = maturity_date(iborIndex, fixing_val_date)
  else
    next_fixing = advance(-Base.Dates.Day(fixingDays), fixing_cal, endDate, iborIndex.convention)
    fixing_end_date = advance(Base.Dates.Day(idx_fixing_days), fixing_cal, next_fixing, iborIndex.convention)
  end

  spanning_time = year_fraction(iborIndex.dc, fixing_val_date, fixing_end_date)

  ## TODO ensure positive (> 0) spanning_time

  return IborCoupon(paymentDate, nominal, startDate, endDate, _fixing_date, fixing_val_date, fixing_end_date, fixingDays, iborIndex, gearing, spread,
                    refPeriodStart, refPeriodEnd, dc, isInArrears, spanning_time)
end


amount(coup::FixedRateCoupon) =
        coup.nominal * (compound_factor(coup.rate, coup.accrualStartDate, coup.accrualEndDate, coup.refPeriodStart, coup.refPeriodEnd) - 1)

amount(coup::IborCoupon, yts::YieldTermStructure) = calc_rate(coup, yts, cap_vol)
date(coup::Coupon) = coup.paymentDate

# legs to build cash flows
abstract Leg <: CashFlows

type FixedRateLeg <: Leg
  coupons::Vector{FixedRateCoupon}
  redemption::SimpleCashFlow

  function FixedRateLeg{B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::B, paymentConvention::C, dc::DC)
    coups = Vector{FixedRateCoupon}(length(schedule.dates) - 1)

    start_date = schedule.dates[1]
    end_date = schedule.dates[2]
    #TODO: setup payment adjustments and the like
    coups[1] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc)

    # build coupons
    count = 2
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]

    while start_date < schedule.dates[end]
      @inbounds coups[count] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc)

      count += 1
      start_date = end_date
      end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    end

    redemption = SimpleCashFlow(faceAmount, end_date)

    new(coups, redemption)
  end
end

type ZeroCouponLeg <: Leg
  redemption::SimpleCashFlow
end

type IborLeg <: Leg
  coupons::Vector{IborCoupon}

  function IborLeg{DC <: DayCount, C <: BusinessDayConvention, I <: Integer, ICP <: IborCouponPricer}(schedule::Schedule, nominal::Float64, iborIndex::IborIndex, paymentDC::DC, paymentAdj::C,
                   fixingDays::Vector{I} = fill(iborIndex.fixingDays, length(schedule.dates) - 1), gearings::Vector{Float64} = ones(length(schedule.dates) - 1),
                   spreads::Vector{Float64} = zeros(length(schedule.dates) - 1), caps::Vector{Float64} = Vector{Float64}(),
                   floors::Vector{Float64} = Vector{Float64}(), isInArrears::Bool = false, isZero::Bool = false,
                   pricer::ICP = BlackIborCouponPricer())
    n = length(schedule.dates) - 1
    coups = Vector{IborCoupon}(n)
    last_payment_date = adjust(schedule.cal, paymentAdj, schedule.dates[end])

    ref_date = _start = end_date = _end = Date()

    for i = 1:n
      ref_start = _start = schedule.dates[i]
      ref_end = _end = schedule.dates[i + 1]
      payment_date = isZero ? last_payment_date : adjust(schedule.cal, paymentAdj, _end)

      # Just Floating Rate Coupons right now
      coups[i] = IborCoupon(payment_date, nominal, _start, _end, fixingDays[i], iborIndex, gearings[i], spreads[i], ref_start, ref_end,
                            paymentDC, isInArrears, pricer)
    end

    new(coups)
  end
end

## Function wrapper for solvers ##
type IRRFinder{L <: Leg, DC <: DayCount, C <: CompoundingType, F <: Frequency}
  leg::L
  npv::Float64
  dc::DC
  comp::C
  freq::F
  includeSettlementDateFlows::Bool
  settlementDate::Date
  npvDate::Date
end

## this function can pass itself or its derivative ##
function operator(finder::IRRFinder)
  function _inner(y::Float64)
    yld = InterestRate(y, finder.dc, finder.comp, finder.freq)
    NPV = npv(finder.leg, yld, finder.includeSettlementDateFlows, finder.settlementDate, finder.npvDate)
    return finder.npv - NPV
  end

  # derivative
  function _inner(y::Float64, ::Derivative)
    yld = InterestRate(y, finder.dc, finder.comp, finder.freq)
    return duration(ModifiedDuration(), finder.leg, yld, finder.dc, finder.includeSettlementDateFlows, finder.settlementDate, finder.npvDate)
  end

  return _inner
end

## Pricer Methods ##
function initialize!(pricer::IborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure)
  # stuff
  if !pricer.initialized
    payment_date = date(coupon)
    if payment_date > yts.referenceDate
      pricer.discount = discount(yts, payment_date)
    else
      pricer.discount = 1.0
    end
    pricer.spreadLegValue = coupon.spread * accrual_period(coupon) * pricer.discount

    pricer.initialized = true
  end

  return pricer
end

function adjusted_fixing(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure,
                        fixing::Float64 = index_fixing(coupon, yts))
  # stuff
  if !coupon.isInArrears
    return fixing
  end

  d1 = coupon.fixingDate
  ref_date = cap_vol.referenceDate

  if d1 <= ref_date
    return fixing
  end

  d2 = value_date(coupon.iborIndex, d1)
  d3 = maturity_date(coupon.iborIndex, d2)
  tau = year_fraction(coupon.iborIndex.dc, d2, d3)
  varience = black_varience(cap_vol, d1, fixing)
  adjustment = fixing * fixing * varience * tau / (1.0 + fixing * tau)
  return fixing + adjustment
end

function swaplet_price(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure)
  swapl_price = adjusted_fixing(pricer, coupon, yts, cap_vol) * accrual_period(coupon) * pricer.discount
  return coupon.gearing * swapl_price + pricer.spreadLegValue
end

swaplet_rate(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure) =
  swaplet_price(pricer, coupon, yts) / (accrual_period(coupon) * pricer.discount)


## Floating Rate Coupon Methods ##
function calc_rate(coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure)
  # first initialize pricer
  initialize!(coupon.pricer, coupon, yts)

  # then get swaplet rate
  return swaplet_rate(coupon.pricer, coupon, yts, cap_vol)
end

function index_fixing(coupon::IborCoupon, yts::YieldTermStructure)
  today = settings.evaluation_date

  if coupon.fixingDate > today
    return forecast_fixing(coupon.iborIndex, yts, coupon.fixingValueDate, coupon.fixingEndDate, coupon.spanningTime)
  end

  error("Fixing date on or before eval date")
end

## NPV METHODS ##
function npv{Y <: YieldTermStructure}(leg::FixedRateLeg, yts::Y, settlement_date::Date, npv_date::Date)
  # stuff
  totalNPV = 0.0
  if length(leg.coupons) == 0
    return totalNPV
  end

  for i in leg
    if has_occurred(i, settlement_date)
      continue
    end

    @inbounds totalNPV += amount(i) * discount(yts, date(i))
  end

  # redemption - not needed with new iterator
  # totalNPV += amount(leg.redemption) * discount(yts, leg.redemption.date)
  return totalNPV / discount(yts, npv_date)
end

function npv(leg::FixedRateLeg, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)
  if length(leg.coupons) == 0
    return 0.0
  end

  totalNPV = 0.0
  discount = 1.0
  last_date = npv_date

  for cp in leg
    if has_occurred(cp, settlement_date)
      continue
    end
    coupon_date = date(cp)
    amount_ = amount(cp)

    if isa(cp, FixedRateCoupon)
      ref_start_date = cp.refPeriodStart
      ref_end_date = cp.refPeriodEnd
    else
      if last_date == npv_date
        ref_start_date = coupon_date - Dates.Year(1)
      else
        ref_start_date = last_date
      end
      ref_end_date = coupon_date
    end

    b = discount_factor(y, last_date, coupon_date, ref_start_date, ref_end_date)

    discount *= b
    last_date = coupon_date

    totalNPV += amount_ * discount
  end

  # redemption - not needed with iterator
  #redempt_date = leg.redemption.date
  # amount = amount(leg.redemption)

  # b = discount_factor(leg.redemption, last_date, redempt_date, last_date, redempt_date)
  # discount *= b
  # totalNPV += amount * discount

  # now we return total npv
  return totalNPV
end

function npv{Y <: YieldTermStructure}(leg::ZeroCouponLeg, yts::Y, settlement_date::Date, npv_date::Date)
  if amount(leg.redemption) == 0.0
    return 0.0
  end

  totalNPV = amount(leg.redemption) * discount(yts, date(leg.redemption))
  return totalNPV / discount(yts, npv_date)
end

## Duration Calculations ##
modified_duration_calc{F <: Frequency}(::SimpleCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::F) = c * B * B * t
modified_duration_calc{F <: Frequency}(::CompoundedCompounding, c::Float64, B::Float64, t::Float64, r::Float64, N::F) = c * t * B / (1 + r / QuantJulia.Time.value(N))
modified_duration_calc{F <: Frequency}(::ContinuousCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::F) = c * B * t
modified_duration_calc{F <: Frequency}(::SimpleThenCompounded, c::Float64, B::Float64, t::Float64, r::Float64, N::F) =
  t <= 1.0 / N ? modified_duration_calc(Simple(), c, B, t, r, N) : modified_duration_calc(CompoundedCompounding(), c, B, t, r, N)

function duration{L <: Leg, DC <: DayCount}(::ModifiedDuration, leg::L, y::InterestRate, dc::DC, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())
  if length(leg.coupons) == 0 # TODO make this applicable to redemption too
    return 0.0
  end

  if npv_date == Date()
    npv_date = settlement_date
  end

  P = 0.0
  t = 0.0
  dPdy = 0.0
  r = y.rate
  N = y.freq
  last_date = npv_date

  for cp in leg
    if has_occurred(cp, settlement_date)
      continue
    end

    c = amount(cp)
    coupon_date = date(cp)
    if isa(cp, FixedRateCoupon)
      ref_start_date = cp.refPeriodStart
      ref_end_date = cp.refPeriodEnd
    else
      if last_date == npv_date
        ref_start_date = coupon_date - Dates.Year(1)
      else
        ref_start_date = last_date
      end
      ref_end_date = coupon_date
    end

    t += year_fraction(dc, last_date, coupon_date, ref_start_date, ref_end_date)

    B = discount_factor(y, t)
    P += c * B

    dPdy -= modified_duration_calc(y.comp, c, B, t, r, N)

    last_date = coupon_date
  end

  if P == 0.0
    return 0.0 # no cashflows
  end

  return -dPdy / P # reverse derivative sign
end

# functions for sorting and finding
sort_cashflow{C <: Coupon}(cf::C) = cf.paymentDate
prev_cf{C <: Coupon}(cf::C, d::Date) = d > cf.paymentDate
next_cf{C <: Coupon}(cf::C, d::Date) = d < cf.paymentDate

function previous_cashflow_date(cf::FixedRateLeg, settlement_date::Date)
  # right now we can assume cashflows are sorted by date because of schedule
  prev_cashflow_idx = findprev(prev_cf, cf.coupons, length(cf.coupons), settlement_date)

  return prev_cashflow_idx == 0 ? 0 : cf.coupons[prev_cashflow_idx].paymentDate
end

function accrual_days{C <: CashFlows, DC <: DayCount}(cf::C, dc::DC, settlement_date::Date)
  last_payment = previous_cashflow_date(cf, settlement_date)

  if last_payment == 0
    return 0.0
  else
    return day_count(dc, last_payment, settlement_date)
  end
end

function next_cashflow{L <: Leg}(cf::L, settlement_date::Date)
  if settlement_date > cf.coupons[end].paymentDate
    return length(cf.coupons)
  end

  return findnext(next_cf, cf.coupons, 1, settlement_date)
end

function accrued_amount(cf::FixedRateLeg, settlement_date::Date, include_settlement_cf::Bool = false)
  next_cf_idx = next_cashflow(cf, settlement_date)
  if cf.coupons[next_cf_idx] == length(cf.coupons)
    return 0.0
  end

  next_cf = cf.coupons[next_cf_idx]
  paymentDate = next_cf.paymentDate
  result = 0.0
  i = next_cf_idx
  while i < length(cf.coupons) && next_cf.paymentDate == paymentDate
    result += accrued_amount(next_cf, settlement_date)
    i += 1
    next_cf = cf.coupons[i]
  end

  return result
end

accrued_amount(cf::ZeroCouponLeg, ::Date, ::Bool= false) = 0.0

function accrued_amount(coup::FixedRateCoupon, settlement_date::Date)
  if settlement_date <= coup.accrualStartDate || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal *
      (compound_factor(coup.rate, coup.accrualStartDate, min(settlement_date, coup.accrualEndDate), coup.refPeriodStart, coup.refPeriodEnd) - 1.0)
end

accrual_period{C <: Coupon}(coup::C) = year_fraction(coup.dc, coup.accrualStartDate, coup.accrualEndDate, coup.refPeriodStart, coup.refPeriodEnd)

function has_occurred{C <: CashFlow}(cf::C, ref_date::Date, include_settlement_cf::Bool = true)
  # will need to expand this
  if ref_date < date(cf) || (ref_date == date(cf) && include_settlement_cf)
    return false
  else
    return true
  end
end

# Yield Calculations ##
function yield{L <: Leg, DC <: DayCount, C <: CompoundingType, F <: Frequency, I <: Integer}(leg::L, npv::Float64, dc::DC, compounding::C, freq::F, include_settlement_cf::Bool, settlement_date::Date,
              npv_date::Date, accuracy::Float64, max_iter::I, guess::Float64)
  solver = NewtonSolver(max_iter)
  obj_fun = IRRFinder(leg, npv, dc, compounding, freq, include_settlement_cf, settlement_date, npv_date)

  return solve(solver, operator(obj_fun), accuracy, guess, guess / 10.0)
end

## ITERATORS ##
# this is to iterate through cash flows and redemption
Base.start(f::FixedRateLeg) = 1
function Base.next(f::FixedRateLeg, state)
  if state > length(f.coupons)
    f.redemption, state + 1
  else
    f.coupons[state], state + 1
  end
end

Base.done(f::FixedRateLeg, state) = length(f.coupons) + 1 < state

# end
