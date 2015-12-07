# building fixed rate bonds
include("github_repos/QuantJulia.jl/src/QuantJulia.jl")
using QuantJulia

function build_bonds(bond_mats, bond_rates, tenor, conv, rule, calendar, dc, freq, issue_date)
  bonds = Vector{FixedRateBond}(10)
  pricing_engine = DiscountingBondEngine()

  for i =1:10
    term_date = bond_mats[i]
    rate = bond_rates[i] / 100
    sched = QuantJulia.Time.Schedule(issue_date, term_date, tenor, conv, conv, rule, true)
    bond = FixedRateBond(0, 100.0, sched, rate, dc, conv, 100.0, issue_date, calendar, pricing_engine)
    bonds[i] = bond
  end

  return bonds
end

function get_npvs(bonds, issue_date, calendar, dc, freq)
  pricing_engine = DiscountingBondEngine()
  rate_quote = Quote(0.05)
  compounding = CompoundedCompounding()
  yts = FlatForwardTermStructure(0, issue_date, calendar, rate_quote, dc, compounding, freq)
  npvs = zeros(length(bonds))
  for i=1:length(bonds)
    npvs[i] = calculate(pricing_engine, yts, bonds[i])
  end

  return npvs
end

function main()
  issue_date = Date(2015, 7, 1)
  bond_mats = [issue_date + Dates.Month(i) for i in range(6, 6, 10)]
  bond_rates = [5.75, 6.0, 6.25, 6.5, 6.75, 6.80, 7.00, 7.1, 7.15, 7.2]
  set_eval_date!(settings, issue_date)

  freq = QuantJulia.Time.Semiannual()
  tenor = QuantJulia.Time.TenorPeriod(freq)
  conv = QuantJulia.Time.Unadjusted()
  rule = QuantJulia.Time.DateGenerationBackwards()
  calendar = QuantJulia.Time.USGovernmentBondCalendar()
  dc = QuantJulia.Time.BondThirty360()
  bonds = build_bonds(bond_mats, bond_rates, tenor, conv, rule, calendar, dc, freq, issue_date)
  npvs = get_npvs(bonds, issue_date, calendar, dc, freq)
end

interp = QuantJulia.Math.LogInterpolation()
trait = Discount()
bootstrap = IterativeBootstrap()

yts = PiecewiseYieldCurve(issue_date, bonds, dc, interp, trait, 0.00000000001, bootstrap)

solver = QuantJulia.Math.BrentSolver()
solver2 = QuantJulia.Math.FiniteDifferenceNewtonSafe()
calculate!(IterativeBootstrap(), yts, solver2, solver)
