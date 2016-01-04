using QuantJulia.Time

function implied_quote{T <: TermStructure}(depo::DepositRate, ts::T)
  return fixing(depo.iborIndex, ts, depo.fixingDate, true)
end
