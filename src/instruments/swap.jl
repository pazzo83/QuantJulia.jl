using QuantJulia.Time

type Payer <: SwapType end
type Receiver <: SwapType end

type VanillaSwap <: Swap
  swapT::SwapType
  nominal::Float64
  fixedSchedule::Schedule
  fixedRate::Float64
  fixedDayCount::DayCount
  iborIndex::IborIndex
  spread::Float64
  floatSchedule::Schedule
  floatDayCount::DayCount
  paymentConvention::BusinessDayConvention
  legs::Vector{Leg}
  payer::Vector{Float64}
end

# Constructors
function VanillaSwap(rate::Float64, tenor::Base.Dates.Period, cal::BusinessCalendar, fixedFrequency::Frequency, fixedConvention::BusinessDayConvention,
                    fixedDayCount::DayCount, iborIndex::IborIndex, spread::Float64, fwdStart::Base.Dates.Period,
                    settlementDays::Integer = iborIndex.fixingDays, nominal::Float64 = 1.0, swapT::SwapType = Payer())
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

  return VanillaSwap(swapT, nominal, fixed_schedule, fixed_rate, fixedDayCount, iborIndex, spread, float_schedule, floatDayCount, fixedConvention,
                    legs, payer)
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
