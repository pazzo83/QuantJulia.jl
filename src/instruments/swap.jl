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
  rate::Quote
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
function VanillaSwap{PrT <: Dates.Period, C <: BusinessCalendar, F <: Frequency, B <: BusinessDayConvention, DC <: DayCount, PrS <: Dates.Period, P <: PricingEngine, I <: Integer, ST <: SwapType}(rate::Float64,
                    tenor::PrT, cal::C, fixedFrequency::F, fixedConvention::B, fixedDayCount::DC, iborIndex::IborIndex, spread::Float64, fwdStart::PrS,
                    pricingEngine::P = DiscountingSwapEngine(), settlementDays::I = iborIndex.fixingDays, nominal::Float64 = 1.0, swapT::ST = Payer())
  # do stuff
  fixedCal = cal
  floatingCal = cal
  floatTenor = iborIndex.tenor
  fixedTenor = QuantJulia.Time.TenorPeriod(fixedFrequency)
  fixedTermConvention = fixedConvention
  floatConvention = iborIndex.convention
  floatTermConvention = iborIndex.convention
  fixedRule = DateGenerationBackwards()
  floatRule = DateGenerationBackwards()
  floatDayCount = iborIndex.dc
  fixed_rate = 0.0

  ref_date = adjust(floatingCal, floatConvention, settings.evaluation_date)
  spot_date = advance(Base.Dates.Day(settlementDays), floatingCal, ref_date, floatConvention)
  start_date = adjust(floatingCal, floatConvention, spot_date + fwdStart)
  ## TODO Float end of month (defaults to false)
  end_date = start_date + tenor

  # build schedules
  fixed_schedule = Schedule(start_date, end_date, fixedTenor, fixedConvention, fixedTermConvention, fixedRule, false, fixedCal)
  float_schedule = Schedule(start_date, end_date, floatTenor, floatConvention, floatTermConvention, floatRule, false, floatingCal)
  # build swap cashflows
  legs = Vector{Leg}(2)
  # first leg is fixed
  legs[1] = FixedRateLeg(fixed_schedule, nominal, fixed_rate, fixedCal, floatConvention, fixedDayCount)
  # second leg is floating
  legs[2] = IborLeg(float_schedule, nominal, iborIndex, floatDayCount, floatConvention)

  payer = _build_payer(swapT)

  results = SwapResults(2)

  return VanillaSwap(Quote(rate), swapT, nominal, fixed_schedule, fixed_rate, fixedDayCount, iborIndex, spread, float_schedule, floatDayCount, fixedConvention,
                    legs, payer, pricingEngine, results, false)
end


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

function maturity_date{S <: Swap}(swap::S)
  d = maturity_date(swap.legs[1])
  for i = 2:length(swap.legs)
    d = max(d, maturity_date(swap.legs[i]))
  end

  return d
end
