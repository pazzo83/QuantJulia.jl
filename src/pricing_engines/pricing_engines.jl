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

type DiscountingSwapEngine{Y <: YieldTermStructure} <: PricingEngine{Y}
  yts::Y

  function call(::Type{DiscountingSwapEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingSwapEngine}, yts::Y)
    new{Y}(yts)
  end
end

function calculate!{B <: Bond}(pe::DiscountingBondEngine, bond::B, recalculate::Bool=false)
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

function calculate!{S <: Swap}(pe::DiscountingSwapEngine, swap::S, recalculate::Bool = false)
  # stuff
  if !swap.calculated || recalculate
    swap.results.value = 0.0
    yts = pe.yts

    ref_date = yts.referenceDate

    swap.results.npvDateDiscount = discount(yts, ref_date)

    for (i, leg) in enumerate(swap.legs)
      swap.results.legNPV[i], swap.results.legBPS[i] = npvbps(leg, yts, ref_date, ref_date)
      swap.results.legNPV[i] *= swap.payer[i]
      swap.results.legBPS[i] *= swap.payer[i]

      d1 = leg.coupons[1].accrualStartDate
      if d1 >= ref_date
        swap.results.startDiscounts[i] = discount(yts, d1)
      end

      d2 = leg.coupons[end].accrualEndDate
      if (d2 >= ref_date)
        swap.results.endDiscounts[i] = discount(yts, d2)
      end

      swap.results.value += swap.results.legNPV[i]
    end
  end

  return swap
end
# end
