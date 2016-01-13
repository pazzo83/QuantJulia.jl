type SettlementPhysical <: SettlementType end
type SettlementCash <: SettlementType end

type Swaption{E <: Exercise, S <: SettlementType} <: Option
  swap::VanillaSwap
  exercise::E
  delivery::S
end

Swaption{E <: Exercise}(swap::VanillaSwap, exercise::E) = Swaption(swap, exercise, SettlementPhysical())
