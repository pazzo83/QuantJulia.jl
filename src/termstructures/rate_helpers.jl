using QuantJulia.Time

function implied_quote(depo::DepositRate, ts::TermStructure)
  return fixing(depo.iborIndex, ts, depo.fixingDate, true)
end
