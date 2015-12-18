# building fixed rate bonds
include("src/QuantJulia.jl")
using QuantJulia

function build_bonds(bond_mats, bond_rates, tenor, conv, rule, calendar, dc, freq, issue_date)
  bonds = Vector{FixedRateBond}(length(bond_mats))
  pricing_engine = DiscountingBondEngine()

  for i =1:length(bond_mats)
    term_date = bond_mats[i]
    # rate = bond_rates[i] / 100
    rate = bond_rates[i]
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

function setup()
  today = now()
  issue_date = Date(Dates.Year(today), Dates.Month(today), Dates.Day(today))
  bond_mats = [issue_date + Dates.Year(i) for i in range(2, 2, 15)]
  # bond_rates = [5.75, 6.0, 6.25, 6.5, 6.75, 6.80, 7.00, 7.1, 7.15, 7.2]
  bond_rates = [0.0200, 0.0225, 0.0250, 0.0275, 0.0300, 0.0325, 0.0350, 0.0375, 0.0400, 0.0425, 0.0450, 0.0475, 0.0500, 0.0525, 0.0550]
  set_eval_date!(settings, issue_date)

  freq = QuantJulia.Time.Annual()
  tenor = QuantJulia.Time.TenorPeriod(freq)
  conv = QuantJulia.Time.Unadjusted()
  rule = QuantJulia.Time.DateGenerationBackwards()
  calendar = QuantJulia.Time.USGovernmentBondCalendar()
  dc = QuantJulia.Time.SimpleDayCount()
  bonds = build_bonds(bond_mats, bond_rates, tenor, conv, rule, calendar, dc, freq, issue_date)

  return issue_date, bonds, dc, calendar
end

function piecewise_yld_curve()
  issue_date, bonds, dc, calendar = setup()

  interp = QuantJulia.Math.LogInterpolation()
  trait = Discount()
  bootstrap = IterativeBootstrap()

  yts = PiecewiseYieldCurve(issue_date, bonds, dc, interp, trait, 0.00000000001, bootstrap)

  solver = QuantJulia.Math.BrentSolver()
  solver2 = QuantJulia.Math.FiniteDifferenceNewtonSafe()
  calculate!(IterativeBootstrap(), yts, solver2, solver)

  println(yts.data)

  # npvs = get_npvs(bonds, issue_date, calendar, dc, freq)
end

function fitted_bond_curve()
  issue_date, bonds, dc, calendar = setup()

  esf = ExponentialSplinesFitting(true, length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, esf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  # println(fitted_curve.fittingMethod.minimumCostValue)
  # println(fitted_curve.fittingMethod.numberOfIterations)
  # println(fitted_curve.fittingMethod.guessSolution)

  disc = discount(fitted_curve, 1.0)
  println(disc)
end

# interp = QuantJulia.Math.LogInterpolation()
# trait = Discount()
# bootstrap = IterativeBootstrap()
#
# yts = PiecewiseYieldCurve(issue_date, bonds, dc, interp, trait, 0.00000000001, bootstrap)
#
# solver = QuantJulia.Math.BrentSolver()
# solver2 = QuantJulia.Math.FiniteDifferenceNewtonSafe()
# calculate!(IterativeBootstrap(), yts, solver2, solver)
#
#
# # fitted bonds
# esf = ExponentialSplinesFitting(true, length(bonds))
# fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, esf, 1e-10, 5000, 1.0)
# initialize!(fitted_curve)
# calculate!(fitted_curve)
