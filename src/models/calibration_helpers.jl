using QuantJulia.Time

type CalibrationHelperCommon
  marketValue::Float64

  CalibrationHelperCommon() = new(0.0)
end

type SwaptionHelper{Dm <: Dates.Period, D1 <: Dates.Period, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure} <: CalibrationHelper
  lazyMixin::LazyMixin
  exerciseDate::Date
  endDate::Date
  maturity::Dm
  swapLength::D1
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

  SwaptionHelper(exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::D1, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TenorPeriod, fixedLegDayCount::DC_fix,
                floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64, yts::T) =
                new(LazyMixin(), exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, CalibrationHelperCommon(), yts)
end

SwaptionHelper{Dm <: Dates.Period, D1 <: Dates.Period, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure}(maturity::Dm, swapLength::D1, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TenorPeriod,
              fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64, yts::T) =
              SwaptionHelper{Dm, Dl, DC_fix, DC_float, T}(Date(), Date(), maturity, swapLength, volatility, index, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, yts)

function perform_calculations!(SwaptionHelper::SwaptionHelper)
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

  return swaptionHelper
end

function add_times_to!(swaptionHelper::SwaptionHelper)
  calculate!(swaptionHelper)

  return swaptionHelper
end

function black_price(swaptionHelper::SwaptionHelper, sigma::Float64)
  calculate!(swaptionHelper)
  # stuff
end
