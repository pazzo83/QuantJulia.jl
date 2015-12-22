using QuantJulia.Time

type FixedRateBond <: Bond
  price::Quote
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

function FixedRateBond(price::Quote, settlementDays::Int64, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DayCount, paymentConvention::BusinessDayConvention,
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

  return FixedRateBond(price, settlementDays, faceAmount, schedule, coups, dc, paymentConvention, redemption, calendar, issueDate, schedule.dates[1], maturityDate, false, pricing_engine, 0.0)
end

type ZeroCouponBond <: Bond
  settlementDays::Integer
  calendar::BusinessCalendar
  faceAmount::Float64
  maturityDate::Date
  paymentConvention::BusinessDayConvention
  redemption::Float64
  cashflows::ZeroCouponLeg
  issueDate::Date
  settlementValue::Float64
  calculated::Bool
end

function ZeroCouponBond(settlementDays::Integer, calendar::BusinessCalendar, faceAmount::Float64, maturityDate::Date, paymentConvention::BusinessDayConvention=Following(),
                        redemption::Float64=100.0, issueDate::Date=Date())
  # build redemption CashFlow
  redemption_cf = ZeroCouponLeg(SimpleCashFlow(redemption, maturityDate))
  return ZeroCouponBond(settlementDays, calendar, faceAmount, maturityDate, paymentConvention, redemption, redemption_cf, issueDate, 0.0, false)
end

value(b::FixedRateBond) = b.price.value
get_settlement_date(b::Bond) = b.issueDate

function notional(bond::Bond, d::Date)
  if d > bond.maturityDate
    return 0.0
  else
    return bond.faceAmount
  end
end

# bond functions
accrued_amount(bond::Bond, settlement::Date) = accrued_amount(bond.cashflows, settlement, false) * 100.0 / notional(bond, settlement)

function yield(bond::Bond, clean_price::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10,
              max_iter::Integer = 100, guess::Float64 = 0.05)
  dirty_price = clean_price + accrued_amount(bond, settlement)
  dirty_price /= 100.0 / notional(bond, settlement)

  return yield(bond.cashflows, dirty_price, dc, compounding, freq, false, settlement, settlement, accuracy, max_iter, guess)
end

function duration(bond::Bond, yld::InterestRate, duration_::Duration, dc::DayCount, settlement_date::Date)
  return duration(duration_, bond.cashflows, yld, dc, false, settlement_date)
end

function duration(bond::Bond, yld::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, duration_::Duration, settlement_date::Date)
  y = InterestRate(yld, dc, compounding, freq)
  return duration(bond, y, duration_, dc, settlement_date)
end

function npv(bond::Bond, pe::PricingEngine, yts::YieldTermStructure)
  calculate!(pe, yts, bond)
  return bond.settlementValue
end
