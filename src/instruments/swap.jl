using QuantJulia.Time

type Payer <: SwapType end
type Receiver <: SwapType end

type SwapResults <: Results
  legNPV::Vector{Float64}
  legBPS::Vector{Float64}
  npvDateDiscount::Float64
  startDiscounts::Vector{Float64}
  endDiscounts::Vector{Float64}
  fairRate::Float64
  value::Float64

  function SwapResults{I <: Integer}(n::I)
    legNPV = zeros(n)
    legBPS = zeros(n)
    startDiscounts = zeros(n)
    endDiscounts = zeros(n)

    new(legNPV, legBPS, 0.0, startDiscounts, endDiscounts, -1.0, 0.0)
  end
end

function reset!(sr::SwapResults)
  n = length(sr.legNPV)
  sr.legNPV = zeros(n)
  sr.legBPS = zeros(n)
  sr.npvDateDiscount = 0.0
  sr.startDiscounts = zeros(n)
  sr.endDiscounts = zeros(n)
  sr.value = 0.0
  sr.fairRate = 0.0

  return sr
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
function VanillaSwap{ST <: SwapType, DC_fix <: DayCount, DC_float <: DayCount, B <: BusinessDayConvention, P <: PricingEngine}(swapT::ST, nominal::Float64, fixedSchedule::Schedule, fixedRate::Float64,
                    fixedDayCount::DC_fix, iborIndex::IborIndex, spread::Float64, floatSchedule::Schedule, floatDayCount::DC_float, pricingEngine::P, paymentConvention::B = floatSchedule.convention)
  # build swap cashflows
  legs = Vector{Leg}(2)
  # first leg is fixed
  legs[1] = FixedRateLeg(fixedSchedule, nominal, fixedRate, fixedSchedule.cal, paymentConvention, fixedDayCount; add_redemption=false)
  # second leg is floating
  legs[2] = IborLeg(floatSchedule, nominal, iborIndex, floatDayCount, paymentConvention; add_redemption=false)

  payer = _build_payer(swapT)

  results = SwapResults(2)

  return VanillaSwap(swapT, nominal, fixedSchedule, fixedRate, fixedDayCount, iborIndex, spread, floatSchedule, floatDayCount, paymentConvention, legs, payer, pricingEngine, results, false)
end

# Swap Helper methods
function _build_payer(swapT::Payer)
  x = ones(2)
  x[1] = -1.0
  return x
end

function _build_payer(swapT::Receiver)
  x = ones(2)
  x[2] = -1.0
  return x
end

# Swap methods #
function maturity_date{S <: Swap}(swap::S)
  d = maturity_date(swap.legs[1])
  for i = 2:length(swap.legs)
    d = max(d, maturity_date(swap.legs[i]))
  end

  return d
end

# Calculation method #
function perform_calculations!(swap::VanillaSwap)
  reset!(swap.results) # reset - TODO this will be expanded
  _calculate!(swap.pricingEngine, swap)

  return swap
end

floating_leg_NPV(swap::VanillaSwap) = swap.results.legNPV[2]
floating_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[2]

fixed_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[1]

function fair_rate(swap::VanillaSwap)
  calculate!(swap)

  return swap.results.fairRate
end
