using QuantJulia.Time

function implied_quote{T <: TermStructure}(depo::DepositRate, ts::T)
  return fixing(depo.iborIndex, ts, depo.fixingDate, true)
end

function implied_quote{S <: Swap}(swap::S)
  const basisPoint = 0.0001
  calculate!(swap.pricingEngine, swap, true)

  floatingLegNPV = floating_leg_NPV(swap)
  spread = swap.spread
  spreadNPV = floating_leg_BPS(swap) / basisPoint * spread
  totNPV = -(floatingLegNPV + spreadNPV)

  return totNPV / (fixed_leg_BPS(swap) / basisPoint)
end
