# bond helper functions
function implied_quote(bond::Bond, yts::YieldTermStructure, pe::PricingEngine, recalculate::Bool=false, clean::Bool = true)
  calculate!(pe, yts, bond, recalculate)
  settlement_value = bond.settlementValue
  return clean ? clean_price(bond, settlement_value, settlement_date(bond)) : dirty_price(bond, settlement_value, settlement_date(bond))
end

clean_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = dirty_price(bond, settlement_value, settlement_date) - accrued_amount(bond, settlement_date)

dirty_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = settlement_value * 100.0 / bond.faceAmount # replace with notionals

function settlement_date(bond::Bond, d::Date = Date())
  if d == Date()
    return settings.evaluation_date
  end

  return d + bond.settlementDays
end

function accrued_amount(bond::Bond, settlement_date::Date)
  return accrued_amount(bond.cashflows, settlement_date) * 100.0 / bond.faceAmount
end
