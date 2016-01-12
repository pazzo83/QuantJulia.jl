using QuantJulia.Time

type Payer <: SwapType end
type Receiver <: SwapType end

type SwapResults <: Results
  legNPV::Vector{Float64}
  legBPS::Vector{Float64}
  npvDateDiscount::Float64
  startDiscounts::Vector{Float64}
  endDiscounts::Vector{Float64}
  value::Float64

  function SwapResults{I <: Integer}(n::I)
    legNPV = zeros(n)
    legBPS = zeros(n)
    startDiscounts = zeros(n)
    endDiscounts = zeros(n)

    new(legNPV, legBPS, 0.0, startDiscounts, endDiscounts, 0.0)
  end
end

type VanillaSwap{ST <: SwapType, DC_fix <: DayCount, DC_float <: DayCount, B <: BusinessDayConvention, L <: Leg, P <: PricingEngine} <: Swap
  swapT::ST
  nominal::Float64
  fixedSchedule::Schedule
  fixedRate::Float64
  fixedDayCount::DC_fix
  iborIndex::IborIndex
  spread::Float64
  floatSchedule::Schedule
  floatDayCount::DC_float
  paymentConvention::B
  legs::Vector{L}
  payer::Vector{Float64}
  pricingEngine::P
  results::SwapResults
  calculated::Bool
end

# Constructors

function maturity_date{S <: Swap}(swap::S)
  d = maturity_date(swap.legs[1])
  for i = 2:length(swap.legs)
    d = max(d, maturity_date(swap.legs[i]))
  end

  return d
end

floating_leg_NPV(swap::VanillaSwap) = swap.results.legNPV[2]
floating_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[2]

fixed_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[1]
