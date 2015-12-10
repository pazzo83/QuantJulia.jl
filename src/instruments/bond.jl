using QuantJulia.Time

type FixedRateBond <: Bond
  settlementDays::Int64
  faceAmount::Float64
  schedule::Schedule
  cashflows::FixedRateLeg
  dc::DayCount
  paymentConvention::BusinessDayConvention
  redemption::Float64
  calendar::BusinessCalendar
  issueDate::Date
  startDate::Date
  maturityDate::Date
  calculated::Bool
  pricing_engine::PricingEngine
  settlementValue::Float64
end

function FixedRateBond(settlementDays::Int64, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DayCount, paymentConvention::BusinessDayConvention,
  redemption::Float64, issueDate::Date, calendar::BusinessCalendar, pricing_engine::PricingEngine)
  maturityDate = schedule.dates[end]

  # num_payments = length(schedule.dates)
  # coups = zeros(num_payments + 1)
  #
  # # build coupons
  # for i = 1:num_payments
  #   fixed_leg = FixedRateCoupon(schedule, faceAmount, coup_rate, calendar, paymentConvention)
  #   coups[i] = fixed_leg
  # end

  # add redemption
  # coups[end] = SimpleCashFlow(redemption, maturityDate)

  coups = FixedRateLeg(schedule, faceAmount, coup_rate, calendar, paymentConvention, dc)

  return FixedRateBond(settlementDays, faceAmount, schedule, coups, dc, paymentConvention, redemption, calendar, issueDate, schedule.dates[1], maturityDate, false, pricing_engine, 0.0)
end

value(b::FixedRateBond) = b.faceAmount
get_settlement_date(b::Bond) = b.issueDate

function notional(bond::Bond, d::Date)
  if d > bond.maturityDate
    return 0.0
  else
    return bond.faceAmount
  end
end

# bond functions
accrued_amount(bond::Bond, settlement::Date) = accrued_amount(bond.cashflows, false, settlement) * 100.0 / notional(bond, settlement)

function yield(bond::Bond, clean_price::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10,
              max_iter::Integer = 100, guess::Float64 = 0.05)
  dirty_price = clean_price + accrued_amount(bond, settlement)
  dirty_price /= 100.0 / notional(bond, settlement)

  return yield(bond.cashflows, dirty_price, dc, compounding, freq, false, settlement, settlement, accuracy, max_iter, guess)
end

function duration(bond::Bond, yld::InterestRate, duration::Duration, settlement_date::Date)
  return duration(duration, bond.cashflows, yld, false, settlement_date)
end

function duration(bond::Bond, yld::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, duration::Duration, settlement_date::Date)
  y = InterestRate(yld, dc, compounding, freq)
  return duration(bond, y, duration, settlement_date)
end
