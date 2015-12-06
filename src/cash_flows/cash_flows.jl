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

# end
