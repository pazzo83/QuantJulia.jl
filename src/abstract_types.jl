# Lazy Object
abstract LazyObject

function calculate!{L <: LazyObject}(lazy::L)
  if !lazy.calculated
    lazy.calculated = true
    perform_calculations!(lazy)
  end

  return lazy
end

function recalculate!{L <: LazyObject}(lazy::L)
  lazy.calculated = false
  calculate!(lazy)

  return lazy
end

# Exercise
abstract Exercise
abstract EarlyExercise <: Exercise

# Instruments
abstract Instrument <: LazyObject
abstract Bond <: Instrument
abstract AbstractRate <: Instrument
abstract Swap <: Instrument
abstract SettlementType
abstract Option <: Instrument
abstract OptionType
abstract SwapType
abstract Results

# Term Structures
abstract TermStructure <: LazyObject
abstract YieldTermStructure <: TermStructure
abstract Curve <: YieldTermStructure
abstract InterpolatedCurve{I, DC, P, T} <: Curve
abstract VolatilityTermStructure <: TermStructure
abstract OptionletVolatilityStructure <: VolatilityTermStructure
abstract SwaptionVolatilityStructure <: VolatilityTermStructure
abstract BootstrapTrait
abstract Bootstrap
abstract FittingMethod
abstract BootstrapHelper <: LazyObject
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

# Models
abstract CalibrationHelper

# Currencies
abstract AbstractCurrency
