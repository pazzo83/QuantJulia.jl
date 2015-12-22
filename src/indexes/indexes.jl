using QuantJulia.Time

type IborIndex <: InterestRateIndex
  familyName::AbstractString
  tenor::Base.Dates.Period
  fixingDays::Integer
  currency::AbstractCurrency
  fixingCalendar::BusinessCalendar
  convention::BusinessDayConvention
  endOfMonth::Bool
  dc::DayCount
end

fixing_date(idx::InterestRateIndex, d::Date) = advance(Base.Dates.Day(-idx.fixingDays), idx.fixingCalendar, d, idx.convention)
maturity_date(idx::IborIndex, d::Date) = advance(idx.tenor, idx.fixingCalendar, d, idx.convention)
value_date(idx::InterestRateIndex, d::Date) = advance(Base.Dates.Day(idx.fixingDays), idx.fixingCalendar, d, idx.convention)

function fixing(idx::InterestRateIndex, ts::TermStructure, fixing_date::Date, forecast_todays_fixing::Bool=true)
  today = settings.evaluation_date
  if fixing_date > today || (fixing_date == today && forecast_todays_fixing)
    return forecast_fixing(idx, ts, fixing_date)
  end

  error("Not yet implemented for older dates than eval date")
end

function forecast_fixing(idx::IborIndex, ts::TermStructure, fixing_date::Date)
  d1 = value_date(idx, fixing_date)
  d2 = maturity_date(idx, d1)
  t = year_fraction(idx.dc, d1, d2)
  return forecast_fixing(idx, ts, d1, d2, t)
end

function forecast_fixing(idx::IborIndex, ts::TermStructure, d1::Date, d2::Date, t::Float64)
  disc1 = discount(ts, d1)
  disc2 = discount(ts, d2)

  return (disc1 / disc2 - 1.0) / t
end
