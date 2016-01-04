# bond helper functions
function implied_quote{B <: Bond}(bond::B, recalculate::Bool=false, clean::Bool = true)
  calculate!(bond.pricingEngine, bond, recalculate)
  settlement_value = bond.settlementValue
  return clean ? clean_price(bond, settlement_value, settlement_date(bond)) : dirty_price(bond, settlement_value, settlement_date(bond))
end

clean_price{B <: Bond}(bond::B, settlement_value::Float64, settlement_date::Date) = dirty_price(bond, settlement_value, settlement_date) - accrued_amount(bond, settlement_date)

dirty_price{B <: Bond}(bond::B, settlement_value::Float64, settlement_date::Date) = settlement_value * 100.0 / bond.faceAmount # replace with notionals

function settlement_date{B <: Bond}(bond::B, d::Date = Date())
  if d == Date()
    d = settings.evaluation_date
  end

  return d + Base.Dates.Day(bond.settlementDays)
end

function accrued_amount{B <: Bond}(bond::B, settlement_date::Date)
  return accrued_amount(bond.cashflows, settlement_date) * 100.0 / bond.faceAmount
end
