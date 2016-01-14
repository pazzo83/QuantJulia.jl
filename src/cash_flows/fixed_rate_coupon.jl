using QuantJulia.Time

type FixedRateCoupon{DC <: DayCount} <: Coupon
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
  accrualStartDate::Date
  accrualEndDate::Date
  refPeriodStart::Date
  refPeriodEnd::Date
  dc::DC
  accrualPeriod::Float64
end

## COUPON METHODS ##
amount(coup::FixedRateCoupon) =
        coup.nominal * (compound_factor(coup.rate, coup.accrualStartDate, coup.accrualEndDate, coup.refPeriodStart, coup.refPeriodEnd) - 1)

calc_rate(coup::FixedRateCoupon) = coup.rate.rate

type FixedRateLeg <: Leg
  coupons::Vector{Union{FixedRateCoupon, SimpleCashFlow}}
  # redemption::SimpleCashFlow

  function FixedRateLeg{B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::B, paymentConvention::C, dc::DC; add_redemption::Bool = true)
    n = add_redemption ? length(schedule.dates) : length(schedule.dates) - 1
    coups = Vector{Union{FixedRateCoupon, SimpleCashFlow}}(n)

    start_date = schedule.dates[1]
    end_date = schedule.dates[2]
    payment_date = adjust(calendar, paymentConvention, end_date)
    #TODO: setup payment adjustments and the like
    coups[1] = FixedRateCoupon{DC}(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

    # build coupons
    count = 2
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    payment_date = adjust(calendar, paymentConvention, end_date)
    while start_date < schedule.dates[end]
      @inbounds coups[count] = FixedRateCoupon{DC}(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

      count += 1
      start_date = end_date
      end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
      payment_date = adjust(calendar, paymentConvention, end_date)
    end

    if add_redemption
      @inbounds coups[end] = SimpleCashFlow(faceAmount, end_date)
    end

    new(coups)
  end
end

# Coupon methods
function accrued_amount(coup::FixedRateCoupon, settlement_date::Date)
  if settlement_date <= coup.accrualStartDate || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal *
      (compound_factor(coup.rate, coup.accrualStartDate, min(settlement_date, coup.accrualEndDate), coup.refPeriodStart, coup.refPeriodEnd) - 1.0)
end
