# cash_flows.jl
# module CF

using QuantJulia.Time

type SimpleCashFlow <: CashFlow
  amount::Float64
  date::Date
end

amount(cf::SimpleCashFlow) = cf.amount

type FixedRateCoupon <: Coupon
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
  accrualStartDate::Date
  accrualEndDate::Date
end

amount(coup::FixedRateCoupon) = coup.nominal * (compound_factor(coup.rate, coup.accrualStartDate, coup.accrualEndDate) - 1)

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
    coups[1] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date)

    # build coupons
    count = 2
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]

    while start_date < schedule.dates[end]
      coups[count] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date)

      count += 1
      start_date = end_date
      end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    end

    redemption = SimpleCashFlow(faceAmount, end_date)

    new(coups, redemption)
  end
end

function npv(leg::Leg, yts::YieldTermStructure, settlement_date::Date, npv_date::Date)
  # stuff
  totalNPV = 0.0
  if length(leg.coupons) == 0
    return totalNPV
  end

  for i in leg.coupons
    #TODO: check has occurred
    if i.paymentDate > settlement_date
      totalNPV += amount(i) * discount(yts, i.paymentDate)
    end
  end

  # redemption
  totalNPV += amount(leg.redemption) * discount(yts, leg.redemption.date)

  return totalNPV / discount(yts, npv_date)
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

function accrued_amount(cf::Leg, settlement_date::Date)
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


# end
