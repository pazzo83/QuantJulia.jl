using QuantJulia.Time

function implied_quote(ts::TermStructure, depo::DepositRate)
  return fixing(depo.iborIndex, ts, depo.fixingDate, true)
end
