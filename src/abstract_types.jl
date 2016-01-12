# Instruments
abstract Instrument
abstract Bond <: Instrument
abstract AbstractRate <: Instrument
abstract Swap <: Instrument
abstract SwapType
abstract Results

# Term Structures
abstract TermStructure
abstract YieldTermStructure <: TermStructure
abstract Curve <: YieldTermStructure
abstract InterpolatedCurve{I, DC, P, T} <: Curve
abstract VolatilityTermStructure <: TermStructure
abstract OptionletVolatilityStructure <: VolatilityTermStructure
abstract BootstrapTrait
abstract Bootstrap
abstract FittingMethod
abstract BootstrapHelper
abstract BondHelper <: BootstrapHelper
abstract RateHelper <: BootstrapHelper

# Pricing Engines
abstract PricingEngine{Y}

# Cash Flows
abstract CashFlows
abstract CashFlow
abstract Coupon <: CashFlow
abstract Duration
abstract CouponPricer
abstract IborCouponPricer <: CouponPricer

# Indexes
abstract InterestRateIndex

# Currencies
abstract AbstractCurrency
