type SettlementPhysical <: SettlementType end
type SettlementCash <: SettlementType end

type SwaptionResults{S <: AbstractString}
  value::Float64
  additionalResults::Dict{S, Float64}
end

SwaptionResults() = SwaptionResults(0.0, Dict("spreadCorrection" => 0.0, "strike" => 0.0, "atmForward" => 0.0, "annuity" => 0.0, "swapLength" => 0.0, "stdDev" => 0.0, "vega" => 0.0))

type Swaption{E <: Exercise, S <: SettlementType, P <: PricingEngine} <: Option
  swap::VanillaSwap
  exercise::E
  delivery::S
  results::SwaptionResults
  pricingEngine::P

  function call{E, S}(::Type(Swaption), swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults)
    new{E, S, PricingEngine}(swap, exercise, delivery, results)
  end

  function call{E, S, P}(::Type(Swaption), swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults, pricingEngine::P)
    new{E, S, P}(swap, exercise, delivery, results, pricingEngine)
  end
end

Swaption{E <: Exercise}(swap::VanillaSwap, exercise::E) = Swaption(swap, exercise, SettlementPhysical(), SwaptionResults())
