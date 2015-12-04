# TermStructures module
module TermStructures

abstract TermStructure
export TermStructure

# term structures
export TermStructure, YieldTermStructure, FlatForwardTermStructure, JumpDate, JumpTime,
  check_range, max_date, time_from_reference, discount, zero_rate, forward_rate, discount_impl
include("yield_term_structure.jl")

# Curves
export InterpolatedCurve, PiecewiseYieldCurve, max_date, discount_impl
include("curves.jl")

# bootstrapping
export BootstrapTrait, Discount, guess, min_value_after, max_value_after,
Bootstrap, IterativeBootstrap, initialize, calculate!, quote_error
include("bootstrap.jl")

# helpers
export implied_quote, clean_price, dirty_price, settlement_date
include("bond_helpers.jl")

end
