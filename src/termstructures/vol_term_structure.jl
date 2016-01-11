using QuantJulia.Time

type NullOptionVolatilityStructure <: OptionletVolatilityStructure end

type ConstantOptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount} <: OptionletVolatilityStructure
  settlementDays::I
  referenceDate::Date
  calendar::B
  bdc::C
  volatility::Float64
  dc::DC
end

function ConstantOptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(settlementDays::I, calendar::B, bdc::C, volatility::Float64, dc::DC)
  today = settings.evaluation_date
  ref_date = advance(Dates.Day(settlementDays), calendar, today, bdc)
  ConstantOptionVolatility(settlementDays, ref_date, calendar, bdc, volatility, dc)
end

function black_varience(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64)
  v = calc_volatility(ovs, option_date, strike)
  t = time_from_reference(ovs, option_date)
  return v * v * t
end

function calc_volatility(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64)
  # TODO checks - see optionletvolatilitystructure.hpp, volatility
  return volatility_impl(ovs, option_date, strike)
end

volatility_impl(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64) =
              volatility_impl(ovs, time_from_reference(ovs, option_date), strike)

volatility_impl(const_opt_vol::ConstantOptionVolatility, ::Float64, ::Float64) = const_opt_vol.volatility
