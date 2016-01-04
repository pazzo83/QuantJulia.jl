using QuantJulia.Time

type IborIndex{S <: AbstractString, I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure} <: InterestRateIndex
  familyName::S
  tenor::TenorPeriod
  fixingDays::I
  currency::AbstractCurrency
  fixingCalendar::B
  convention::C
  endOfMonth::Bool
  dc::DC
  ts::T

  call{S, I, B, C, DC}(::Type{IborIndex}, familyName::S, tenor::TenorPeriod, fixingDays::I, currency::AbstractCurrency, fixingCalendar::B,
                                convention::C, endOfMonth::Bool, dc::DC) =
    new{S, I, B, C, DC, TermStructure}(familyName, tenor, fixingDays, currency, fixingCalendar, convention, endOfMonth, dc)
end

fixing_date{I <: InterestRateIndex}(idx::I, d::Date) = advance(Dates.Day(-idx.fixingDays), idx.fixingCalendar, d, idx.convention)
maturity_date(idx::IborIndex, d::Date) = advance(idx.tenor.period, idx.fixingCalendar, d, idx.convention)
value_date{I <: InterestRateIndex}(idx::I, d::Date) = advance(Dates.Day(idx.fixingDays), idx.fixingCalendar, d, idx.convention)

function fixing{I <: InterestRateIndex, T <: TermStructure}(idx::I, ts::T, fixing_date::Date, forecast_todays_fixing::Bool=true)
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

# types of indexes
function euribor_index(tenor::TenorPeriod)
  return IborIndex("Euribor", tenor, 2, EURCurrency(), QuantJulia.Time.TargetCalendar(), euribor_conv(tenor.period), euribor_eom(tenor.period),
                  QuantJulia.Time.Actual360())
end

euribor_conv(::Union{Base.Dates.Day, Base.Dates.Week}) = QuantJulia.Time.Following()
euribor_conv(::Union{Base.Dates.Month, Base.Dates.Year}) = QuantJulia.Time.ModifiedFollowing()

euribor_eom(::Union{Base.Dates.Day, Base.Dates.Week}) = false
euribor_eom(::Union{Base.Dates.Month, Base.Dates.Year}) = true
