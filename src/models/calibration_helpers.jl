using QuantJulia.Time

type CalibrationHelperCommon
  marketValue::Float64

  CalibrationHelperCommon() = new(0.0)
end

type SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure} <: CalibrationHelper
  lazyMixin::LazyMixin
  exerciseDate::Date
  endDate::Date
  maturity::Dm
  swapLength::Dl
  volatility::Quote
  iborIndex::IborIndex
  fixedLegTenor::TenorPeriod
  fixedLegDayCount::DC_fix
  floatingLegDayCount::DC_float
  strike::Float64
  nominal::Float64
  shift::Float64
  exerciseRate::Float64
  calibCommon::CalibrationHelperCommon
  yts::T
  swaption::Swaption

  SwaptionHelper(exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TenorPeriod, fixedLegDayCount::DC_fix,
                floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64, yts::T) =
                new(LazyMixin(), exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, CalibrationHelperCommon(), yts)
end

SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure}(maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TenorPeriod,
              fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, yts::T, strike::Float64 = -1.0, nominal::Float64 = 1.0, shift::Float64 = 0.0, exerciseRate::Float64 = 0.0) =
              SwaptionHelper{Dm, Dl, DC_fix, DC_float, T}(Date(), Date(), maturity, swapLength, volatility, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, yts)

function perform_calculations!(swaptionHelper::SwaptionHelper)
  calendar = swaptionHelper.iborIndex.fixingCalendar
  fixingDays = swaptionHelper.iborIndex.fixingDays
  convention = swaptionHelper.iborIndex.convention

  exerciseDate = swaptionHelper.exerciseDate == Date() ? advance(swaptionHelper.maturity, calendar, swaptionHelper.yts.referenceDate, convention) : swaptionHelper.exerciseDate

  startDate = advance(Dates.Day(fixingDays), calendar, exerciseDate, convention)

  endDate = swaptionHelper.endDate == Date() ? advance(swaptionHelper.swapLength, calendar, startDate, convention) : swaptionHelper.endDate

  fixedSchedule = QuantJulia.Time.Schedule(startDate, endDate, swaptionHelper.fixedLegTenor, convention, convention, QuantJulia.Time.DateGenerationForwards(), false, calendar)
  floatSchedule = QuantJulia.Time.Schedule(startDate, endDate, swaptionHelper.iborIndex.tenor, convention, convention, QuantJulia.Time.DateGenerationForwards(), false, calendar)

  swapEngine = DiscountingSwapEngine(swaptionHelper.yts)

  swapT = Receiver()

  tempSwap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, 0.0, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)

  forward = fair_rate(tempSwap)

  if swaptionHelper.strike == -1.0
    swaptionHelper.exerciseRate = forward
  else
    swaptionHelper.exerciseRate = strike
    swapT = strike <= forward ? Receiver() : Payer()
  end

  swap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, swaptionHelper.exerciseRate, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)
  exercise = EuropeanExercise(exerciseDate)

  swaptionHelper.swaption = Swaption(swap, exercise)

  # calibration calc
  _calibration_calculate!(swaptionHelper)

  return swaptionHelper
end

function _calibration_calculate!(swaptionHelper::SwaptionHelper)
  swaptionHelper.calibCommon.marketValue = black_price!(swaptionHelper, swaptionHelper.volatility.value)

  return swaptionHelper
end

function add_times_to!(swaptionHelper::SwaptionHelper, times::Vector{Float64})
  calculate!(swaptionHelper)
  discretizedSwap = DiscretizedSwaption(swaptionHelper.swaption, reference_date(swaptionHelper.yts), swaptionHelper.yts.dc)
  times = vcat(times, mandatory_times(discretizedSwap))
  return times
end

function black_price!(swaptionHelper::SwaptionHelper, sigma::Float64)
  calculate!(swaptionHelper)
  # stuff
  black = BlackSwaptionEngine(swaptionHelper.yts, Quote(sigma), Actual365(), swaptionHelper.shift)
  swaptionHelper.swaption.pricingEngine = black
  value = npv(swaptionHelper.swaption)
  return value
end