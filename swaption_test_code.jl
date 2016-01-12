include("src/QuantJulia.jl")
using QuantJulia

function generate_flatforward_ts()
  cal = QuantJulia.Time.TargetCalendar()
  settlement_date = Date(2002, 2, 19)
  todays_date = Date(2002, 2, 15)
  set_eval_date!(settings, todays_date)

  flat_rate = Quote(0.04875825)

  ffts = FlatForwardTermStructure(settlement_date, cal, flat_rate, QuantJulia.Time.Actual365())

  return ffts
end
