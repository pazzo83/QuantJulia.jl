# Day Count (adapted from Ito.jl and InterestRates.jl)

abstract DayCount

type Actual360 <:DayCount ; end
type Actual365 <: DayCount ; end
abstract Thirty360 <:DayCount

type BondThirty360 <: Thirty360; end
type EuroBondThirty360 <: Thirty360; end
type ItalianThirty360 <: Thirty360; end

typealias USAThirty360 BondThirty360
typealias EuroThirty360 EuroBondThirty360

abstract ActualActual <:DayCount

type ISMAActualActual <: ActualActual; end
type ISDAActualActual <: ActualActual; end
type AFBActualActual <: ActualActual; end

type SimpleDayCount <: DayCount end

# Day Counting
# default day count method
day_count(c::DayCount, d_start::Date, d_end::Date) = Int(d_end - d_start)

# days per year
days_per_year(::Union{Actual360, Thirty360}) = 360
days_per_year(::Actual365) = 365

# year fractions
# default
year_fraction(c::SimpleDayCount, d_start::Date, d_end::Date) = year_fraction(c, d_start, d_end, Date(), Date())

year_fraction(c::DayCount, d_start::Date, d_end::Date) = day_count(c, d_start, d_end) / days_per_year(c)

# add'l methods
year_fraction(c::Union{Actual360, Thirty360, Actual365}, d_start::Date, d_end::Date, ::Date, ::Date) = year_fraction(c, d_start, d_end)



function year_fraction(::SimpleDayCount, d_start::Date, d_end::Date, ::Date, ::Date)
  dm_start = Dates.Day(d_start)
  dm_end = Dates.Day(d_end)

  if dm_start == dm_end || (dm_start > dm_end && Dates.lastdayofmonth(d_end) == d_end) || (dm_start < dm_end && Dates.lastdayofmonth(d_start) == d_start)
    return Int(Dates.Year(d_end) - Dates.Year(d_start)) +  (Int(Dates.Month(d_end)) - Int(Dates.Month(d_start))) / 12.0
  else
    return year_fraction(BondThirty360(), d_start, d_end)
  end
end
