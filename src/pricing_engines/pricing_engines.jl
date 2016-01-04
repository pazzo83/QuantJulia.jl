# Pricing Engines module
# module PricingEngines

type DiscountingBondEngine{Y <: YieldTermStructure} <: PricingEngine{Y}
  yts::Y

  function call(::Type{DiscountingBondEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingBondEngine}, yts::Y)
    new{Y}(yts)
  end
end

function calculate!{E <: DiscountingBondEngine, B <: Bond}(pe::E, bond::B, recalculate::Bool=false)
  if bond.calculated && !recalculate
    return bond
  end
  yts = pe.yts
  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  bond.calculated = true
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

# end
