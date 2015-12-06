# Pricing Engines module
# module PricingEngines

type DiscountingBondEngine <: PricingEngine end

function calculate(::DiscountingBondEngine, yts::YieldTermStructure, bond::Bond)
  valuation_date = yts.reference_date
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)

  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return value
end

# end
