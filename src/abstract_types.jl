# Instruments
abstract Instrument
abstract Bond <: Instrument

# Term Structures
abstract TermStructure
abstract YieldTermStructure <: TermStructure
abstract InterpolatedCurve{I, T, B} <: YieldTermStructure
abstract BootstrapTrait
abstract Bootstrap

# Pricing Engines
abstract PricingEngine

# Cash Flows
abstract CashFlows
abstract CashFlow
abstract Coupon <: CashFlow
