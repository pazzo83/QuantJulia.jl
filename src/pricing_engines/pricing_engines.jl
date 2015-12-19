# Pricing Engines module
# module PricingEngines

type DiscountingBondEngine <: PricingEngine end

function calculate!(::DiscountingBondEngine, yts::YieldTermStructure, bond::Bond, recalculate::Bool=false)
  if bond.calculated && !recalculate
    return bond
  end

  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  bond.calculated = true
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

# end
