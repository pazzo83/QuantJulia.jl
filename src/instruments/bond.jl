using QuantJulia.Time

type FixedRateBond{I <: Integer, DC <: DayCount, B <: BusinessDayConvention, C <: BusinessCalendar, P <: PricingEngine} <: Bond
  price::Quote
  settlementDays::I
  faceAmount::Float64
  schedule::Schedule
  cashflows::FixedRateLeg
  dc::DC
  paymentConvention::B
  redemption::Float64
  calendar::C
  issueDate::Date
  startDate::Date
  maturityDate::Date
  calculated::Bool
  pricingEngine::P
  settlementValue::Float64
end

function FixedRateBond{I <: Integer, DC <: DayCount, B <: BusinessDayConvention, C <: BusinessCalendar, P <: PricingEngine}(price::Quote, settlementDays::I, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DC, paymentConvention::B,
  redemption::Float64, issueDate::Date, calendar::C, pricing_engine::P)
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

  coups = FixedRateLeg{FixedRateCoupon, SimpleCashFlow}(schedule, faceAmount, coup_rate, calendar, paymentConvention, dc)

  return FixedRateBond(price, settlementDays, faceAmount, schedule, coups, dc, paymentConvention, redemption, calendar, issueDate, schedule.dates[1], maturityDate, false, pricing_engine, 0.0)
end

type FloatingRateBond{I <: Integer, DC <: DayCount, B <: BusinessDayConvention, P <: PricingEngine} <: Bond
  settlementDays::I
  faceAmount::Float64
  schedule::Schedule
  cashflows::IborLeg
  iborIndex::IborIndex
  dc::DC
  convention::B
  fixingDays::I
  gearings::Vector{Float64}
  spreads::Vector{Float64}
  caps::Vector{Float64}
  floors::Vector{Float64}
  inArrears::Bool
  redemption::Float64
  issueDate::Date
  maturityDate::Date
  calculated::Bool
  pricingEngine::P
  settlementValue::Float64
end

function FloatingRateBond{I <: Integer, DC <: DayCount, B <: BusinessDayConvention, P <: PricingEngine}(settlementDays::I, faceAmount::Float64, schedule::Schedule, iborIndex::IborIndex, dc::DC,
                          convention::B, fixingDays::I, gearings::Vector{Float64}, spreads::Vector{Float64},
                          caps::Vector{Float64}, floors::Vector{Float64}, inArrears::Bool, redemption::Float64, issueDate::Date,
                          pricingEngine::P)
  maturityDate = schedule.dates[end]
  fixingDaysVect = fill(fixingDays, length(schedule.dates) - 1)

  coups = IborLeg(schedule, faceAmount, iborIndex, dc, convention, fixingDaysVect, gearings, spreads, caps, floors, inArrears)
  return FloatingRateBond(settlementDays, faceAmount, schedule, coups, iborIndex, dc, convention, fixingDays, gearings, spreads, caps, floors,
                          inArrears, redemption, issueDate, maturityDate, false, pricingEngine, 0.0)
end

type ZeroCouponBond{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention} <: Bond
  settlementDays::I
  calendar::B
  faceAmount::Float64
  maturityDate::Date
  paymentConvention::C
  redemption::Float64
  cashflows::ZeroCouponLeg
  issueDate::Date
  settlementValue::Float64
  calculated::Bool
end

function ZeroCouponBond{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention}(settlementDays::I, calendar::B, faceAmount::Float64, maturityDate::Date, paymentConvention::C=Following(),
                        redemption::Float64=100.0, issueDate::Date=Date())
  # build redemption CashFlow
  redemption_cf = ZeroCouponLeg(SimpleCashFlow(redemption, maturityDate))
  return ZeroCouponBond(settlementDays, calendar, faceAmount, maturityDate, paymentConvention, redemption, redemption_cf, issueDate, 0.0, false)
end

value(b::FixedRateBond) = b.price.value
get_settlement_date{B <: Bond}(b::B) = b.issueDate

function notional{B <: Bond}(bond::B, d::Date)
  if d > bond.maturityDate
    return 0.0
  else
    return bond.faceAmount
  end
end

# bond functions
accrued_amount{B <: Bond}(bond::B, settlement::Date) = accrued_amount(bond.cashflows, settlement, false) * 100.0 / notional(bond, settlement)

maturity_date{B <: Bond}(bond::B) = maturity_date(bond.cashflows)

function yield{B <: Bond, DC <: DayCount, C <: CompoundingType, F <: Frequency, I <: Integer}(bond::B, clean_price::Float64, dc::DC, compounding::C, freq::F, settlement::Date, accuracy::Float64 = 1.0e-10,
              max_iter::I = 100, guess::Float64 = 0.05)
  dirty_price = clean_price + accrued_amount(bond, settlement)
  dirty_price /= 100.0 / notional(bond, settlement)

  return yield(bond.cashflows, dirty_price, dc, compounding, freq, false, settlement, settlement, accuracy, max_iter, guess)
end

function duration{B <: Bond, D <: Duration, DC <: DayCount}(bond::B, yld::InterestRate, duration_::D, dc::DC, settlement_date::Date)
  return duration(duration_, bond.cashflows, yld, dc, false, settlement_date)
end

function duration{B <: Bond, DC <: DayCount, C <: CompoundingType, F <: Frequency, D <: Duration}(bond::B, yld::Float64, dc::DC, compounding::C, freq::F, duration_::D, settlement_date::Date)
  y = InterestRate(yld, dc, compounding, freq)
  return duration(bond, y, duration_, dc, settlement_date)
end

function npv{B <: Bond, P <: PricingEngine}(bond::B, pe::P)
  calculate!(pe, bond)
  return bond.settlementValue
end

clean_price{B <: Bond}(bond::B) = clean_price(bond, bond.settlementValue, settlement_date(bond))
dirty_price{B <: Bond}(bond::B) = dirty_price(bond, bond.settlementValue, settlement_date(bond))
