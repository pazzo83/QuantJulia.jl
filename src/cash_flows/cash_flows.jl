# cash_flows.jl
# module CF

using QuantJulia.Time

## TYPES ##
## types of cash flows ##
type SimpleCashFlow <: CashFlow
  amount::Float64
  date::Date
end

amount(cf::SimpleCashFlow) = cf.amount
date(cf::SimpleCashFlow) = cf.date

type FixedRateCoupon <: Coupon
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
  accrualStartDate::Date
  accrualEndDate::Date
  refPeriodStart::Date
  refPeriodEnd::Date
end

amount(coup::FixedRateCoupon) = coup.nominal * (compound_factor(coup.rate, coup.accrualStartDate, coup.accrualEndDate) - 1)
date(coup::FixedRateCoupon) = coup.paymentDate

# legs to build cash flows
abstract Leg <: CashFlows

type FixedRateLeg <: Leg
  coupons::Vector{FixedRateCoupon}
  redemption::SimpleCashFlow

  function FixedRateLeg(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::BusinessCalendar, paymentConvention::BusinessDayConvention, dc::DayCount)
    coups = Vector{FixedRateCoupon}(length(schedule.dates) - 1)

    start_date = schedule.dates[1]
    end_date = schedule.dates[2]
    #TODO: setup payment adjustments and the like
    coups[1] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date)

    # build coupons
    count = 2
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]

    while start_date < schedule.dates[end]
      @inbounds coups[count] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date)

      count += 1
      start_date = end_date
      end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    end

    redemption = SimpleCashFlow(faceAmount, end_date)

    new(coups, redemption)
  end
end

## Function wrapper for solvers ##
type IRRFinder
  leg::Leg
  npv::Float64
  dc::DayCount
  comp::CompoundingType
  freq::Frequency
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

## NPV METHODS ##
function npv(leg::Leg, yts::YieldTermStructure, settlement_date::Date, npv_date::Date)
  # stuff
  totalNPV = 0.0
  if length(leg.coupons) == 0
    return totalNPV
  end

  for i in leg
    #TODO: check has occurred
    @inbounds if date(i) > settlement_date
      @inbounds totalNPV += amount(i) * discount(yts, date(i))
    end
  end

  # redemption - not needed with new iterator
  # totalNPV += amount(leg.redemption) * discount(yts, leg.redemption.date)

  return totalNPV / discount(yts, npv_date)
end

function npv(leg::Leg, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)
  if length(leg.coupons) == 0
    return 0.0
  end

  totalNPV = 0.0
  discount = 1.0
  last_date = npv_date

  for cp in leg
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

## Duration Calculations ##
modified_duration_calc(::SimpleCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::Frequency) = c * B * B * t
modified_duration_calc(::CompoundedCompounding, c::Float64, B::Float64, t::Float64, r::Float64, N::Frequency) = c * t * B / (1 + r / QuantJulia.Time.value(N))
modified_duration_calc(::ContinuousCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::Frequency) = c * B * t
modified_duration_calc(::SimpleThenCompounded, c::Float64, B::Float64, t::Float64, r::Float64, N::Frequency) =
  t <= 1.0 / N ? modified_duration_calc(Simple(), c, B, t, r, N) : modified_duration_calc(CompoundedCompounding(), c, B, t, r, N)

function duration(::ModifiedDuration, leg::Leg, y::InterestRate, dc::DayCount, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())
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
    # TODO check has occurred
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
sort_cashflow(cf::Coupon) = cf.paymentDate
prev_cf(cf::Coupon, d::Date) = d > cf.paymentDate
next_cf(cf::Coupon, d::Date) = d < cf.paymentDate

function previous_cashflow_date(cf::FixedRateLeg, settlement_date::Date)
  # right now we can assume cashflows are sorted by date because of schedule
  prev_cashflow_idx = findprev(prev_cf, cf.coupons, length(cf.coupons), settlement_date)

  return prev_cashflow_idx == 0 ? 0 : cf.coupons[prev_cashflow_idx].paymentDate
end

function accrual_days(cf::CashFlows, dc::DayCount, settlement_date::Date)
  last_payment = previous_cashflow_date(cf, settlement_date)

  if last_payment == 0
    return 0.0
  else
    return day_count(dc, last_payment, settlement_date)
  end
end

function next_cashflow(cf::Leg, settlement_date::Date)
  if settlement_date > cf.coupons[end].paymentDate
    return length(cf.coupons)
  end

  return findnext(next_cf, cf.coupons, 1, settlement_date)
end

function accrued_amount(cf::Leg, settlement_date::Date, include_settlement_cf::Bool = false)
  next_cf_idx = next_cashflow(cf, settlement_date)

  if cf.coupons[next_cf_idx] == length(cf.coupons)
    return 0.0
  end

  result = 0.0

  next_cf = cf.coupons[next_cf_idx]
  paymentDate = next_cf.paymentDate
  i = 0
  while next_cf == paymentDate
    result += accrued_amount(next_cf, settlement_date)
  end

  return result
end

function accrued_amount(coup::FixedRateCoupon, settlement_date::Date)
  if settlement_date <= coup.accrualStartDate || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal * (compound_factor(coup.rate, coup.accrualStartDate, min(settlement_date, coup.accrualEndDate)) - 1.0)
end

function has_occurred(cf::CashFlow, ref_date::Date)
  # will need to expand this
  if ref_date < date(cf)
    return false
  else
    return true
  end
end

# Yield Calculations ##
function yield(leg::Leg, npv::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, include_settlement_cf::Bool, settlement_date::Date,
              npv_date::Date, accuracy::Float64, max_iter::Integer, guess::Float64)
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
