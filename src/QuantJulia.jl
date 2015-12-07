# QuantJulia

module QuantJulia

# functions overridden from based
import Base.findprev, Base.findnext

function findprev(testf::Function, A, start::Integer, val)
  for i = start:-1:1
    testf(A[i], val) && return i
  end
  0
end

function findnext(testf::Function, A, start::Integer, val)
  for i = start:length(A)
    if testf(A[i], val)
      return i
    end
  end
  return 0
end

# Time module
include("time/Time.jl")

# Math module
include("math/Math.jl")

# MAIN MODULE CODE

export
    # abstract_types.jl
    CompoundingType, TermStructure, YieldTermStructure, InterpolatedCurve, BootstrapTrait, Bootstrap,
    CashFlows, CashFlow, Coupon, Instrument, Bond, PricingEngine,

    # quotes/Quotes/jl
    Quote,

    # InterestRates.jl
    ContinuousCompounding, SimpleCompounding, CompoundedCompounding, SimpleThenCompounded,
    InterestRate, discount_factor, compound_factor, equivalent_rate, implied_rate,

    # termstructures/yield_term_structure.jl
    FlatForwardTermStructure, JumpDate, JumpTime,
    check_range, max_date, time_from_reference, discount, zero_rate, forward_rate, discount_impl,

    # termstructures/curve.jl
    PiecewiseYieldCurve, max_date,

    # termstructures/bootstrap.jl
    Discount, guess, min_value_after, max_value_after,
    IterativeBootstrap, initialize, calculate!, quote_error,

    # termstructures/bond_helpers.jl
    implied_quote, clean_price, dirty_price, settlement_date,

    # cash_flows/cash_flows.jl
    SimpleCashFlow, FixedRateCoupon, Leg, FixedRateLeg,

    # instruments/bond.jl
    FixedRateBond,

    # pricing_engines/pricing_engines.jl
    DiscountingBondEngine, calculate

# abstract types
include("abstract_types.jl")

# Quotes ----------------------------
include("quotes/Quotes.jl")

# Interest Rates ---------------------------------
include("InterestRates.jl")

# Term Structures -----------------------------------
# yield term structures
include("termstructures/yield_term_structure.jl")
# Curves
include("termstructures/curve.jl")
# bootstrapping
include("termstructures/bootstrap.jl")
# helpers
include("termstructures/bond_helpers.jl")

# Cash Flows ------------------------------------
include("cash_flows/cash_flows.jl")

# Instruments ------------------------
# bond
include("instruments/bond.jl")

# Pricing Engines ------------------------
include("pricing_engines/pricing_engines.jl")

# # Helpers NOW IN TERM STRUCTURE
# include("helpers/bond_helpers.jl")

type Settings
  evaluation_date::Date
end

settings = Settings(Date())

function set_eval_date!(sett::Settings, d::Date)
  sett.evaluation_date = d
end

export Settings, settings, set_eval_date!

end
