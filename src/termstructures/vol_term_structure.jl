using QuantJulia.Time

type NullOptionVolatilityStructure <: OptionletVolatilityStructure end

type ConstantOptionVolatility <: OptionletVolatilityStructure
  settlementDays::Integer
  referenceDate::Date
  calendar::BusinessCalendar
  bdc::BusinessDayConvention
  volatility::Float64
  dc::DayCount
end

function ConstantOptionVolatility(settlementDays::Integer, calendar::BusinessCalendar, bdc::BusinessDayConvention, volatility::Float64, dc::DayCount)
  today = settings.evaluation_date
  ref_date = advance(Dates.Day(settlementDays), calendar, today)
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

volatility_impl(const_opt_vol::ConstantOptionVolatility, ::Float64, ::Float64) = const_opt_val.volatility
