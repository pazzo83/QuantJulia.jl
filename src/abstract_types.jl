# Instruments
abstract Instrument
abstract Bond <: Instrument
abstract AbstractRate <: Instrument

# Term Structures
abstract TermStructure
abstract YieldTermStructure <: TermStructure
abstract Curve <: YieldTermStructure
abstract InterpolatedCurve{I} <: Curve
abstract BootstrapTrait
abstract Bootstrap
abstract FittingMethod

# Pricing Engines
abstract PricingEngine

# Cash Flows
abstract CashFlows
abstract CashFlow
abstract Coupon <: CashFlow
abstract Duration

# Indexes
abstract InterestRateIndex

# Currencies
abstract AbstractCurrency
