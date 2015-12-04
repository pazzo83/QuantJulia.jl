using QuantJulia.Time, QuantJulia.CF, QuantJulia.TermStructures

abstract Bond <: Instrument

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
end

function FixedRateBond(settlementDays::Int64, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DayCount, paymentConvention::BusinessDayConvention,
  redemption::Float64, issueDate::Date, calendar::BusinessCalendar)
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

  return FixedRateBond(settlementDays, faceAmount, schedule, coups, dc, paymentConvention, redemption, calendar, issueDate, schedule.dates[1], maturityDate, false)
end
