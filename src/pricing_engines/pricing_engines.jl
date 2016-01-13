# Pricing Engines module
# module PricingEngines
using QuantJulia.Time

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

type BlackSwaptionEngine{Y <: YieldTermStructure, S <: SwaptionVolatilityStructure, DC <: DayCount} <: PricingEngine{Y}
  yts::Y
  vol::Quote
  volStructure::S
  dc::DC
  displacement::Float64
end

BlackSwaptionEngine{Y <: YieldTermStructure, DC <: DayCount}(yts::Y, vol::Quote, dc::DC, displacement::Float64 = 0.0) =
                    BlackSwaptionEngine(yts, vol, ConstantSwaptionVolatility(0, QuantJulia.Time.NullCalendar(), QuantJulia.Time.Following(), vol, dc), dc, displacement)


function _calculate!{B <: Bond}(pe::DiscountingBondEngine, bond::B)
  yts = pe.yts
  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

function _calculate!{S <: Swap}(pe::DiscountingSwapEngine, swap::S)
  # stuff
  # println("NEW ONE=============================================================================")
  # if swap.rate.value > 0.0323
  #   error("DUIE")
  # end
  const basisPoint = 0.0001
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

  if swap.results.legBPS[1] != 0.0
    swap.results.fairRate = swap.fixedRate - swap.results.value / (swap.results.legBPS[1] / basisPoint)
  end

  return swap
end
# end
