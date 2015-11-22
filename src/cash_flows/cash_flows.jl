# cash_flows.jl
module CF

using QuantJulia.Time, QuantJulia.InterestRates

export CashFlows, CashFlow, SimpleCashFlow, Coupon, FixedRateCoupon, Leg, FixedRateLeg
# abstract type for all cash flows
abstract CashFlows

# abstract type for a single cashflow leg
abstract CashFlow

type SimpleCashFlow <: CashFlow
  amount::Float64
  date::Date
end

# coupons
abstract Coupon <: CashFlow

type FixedRateCoupon <: Coupon
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
  calendar::BusinessCalendar
  paymentConvention::BusinessDayConvention
end

# legs to build cash flows
abstract Leg <: CashFlows

type FixedRateLeg <: Leg
  coupons::Vector{FixedRateCoupon}
  redemption::SimpleCashFlow

  function FixedRateLeg(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::BusinessCalendar, paymentConvention::BusinessDayConvention, dc::DayCount)
    coups = Vector{FixedRateCoupon}(length(schedule.dates))

    start_date = schedule.dates[1]
    end_date = schedule.dates[2]
    #TODO: setup payment adjustments and the like
    coups[1] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), calendar, paymentConvention)

    for i=2:length(schedule.dates)
      start_date = end_date # we will need this in the future for accrual dates
      end_date = i == length(schedule.dates) ? schedule.dates[end] : schedule.dates[i + 1]

      coups[i] = FixedRateCoupon(end_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), calendar, paymentConvention)
    end

    redemption = SimpleCashFlow(faceAmount, end_date)

    new(coups, redemption)
  end
end

end
