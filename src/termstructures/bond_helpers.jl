# Bond Helpers
type FixedRateBondHelper <: BondHelper
  price::Quote
  bond::FixedRateBond
end
value(b::FixedRateBondHelper) = b.price.value
maturity_date(b::FixedRateBondHelper) = maturity_date(b.bond)

# bond helper functions
function implied_quote{B <: BondHelper}(bond_h::B, recalculate::Bool=false, clean::Bool = true)
  bond = bond_h.bond
  calculate!(bond.pricingEngine, bond, recalculate)
  settlement_value = bond.settlementValue
  return clean ? clean_price(bond, settlement_value, settlement_date(bond)) : dirty_price(bond, settlement_value, settlement_date(bond))
end
