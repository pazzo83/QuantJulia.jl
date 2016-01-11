using QuantJulia.Time

function implied_quote{T <: TermStructure}(depo::DepositRate, ts::T)
  return fixing(depo.iborIndex, ts, depo.fixingDate, true)
end

function implied_quote{S <: Swap}(swap::S)
  const basisPoint = 0.0001
  calculate!(swap.pricingEngine, swap, true)
  #
  # println("Floating Leg NPV ", floating_leg_NPV(swap))
  # println("Floating Leg BPS ", floating_leg_BPS(swap))
  # println("Fixed Leg BPS ", fixed_leg_BPS(swap))
  # println("Swap spread ", swap.spread)

  floatingLegNPV = floating_leg_NPV(swap)
  spread = swap.spread
  spreadNPV = floating_leg_BPS(swap) / basisPoint * spread
  totNPV = -(floatingLegNPV + spreadNPV)

  return totNPV / (fixed_leg_BPS(swap) / basisPoint)
end
