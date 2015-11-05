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

# Day Counting
# default day count method
day_count(c::DayCount, d_start::Date, d_end::Date) = Int(d_end - d_start)

# days per year
days_per_year(::Union{Actual360, Thirty360}) = 360
days_per_year(::Actual365) = 365

# year fractions
# default
year_fraction(c::DayCount, d_start::Date, d_end::Date) = day_count(c, d_start, d_end) / days_per_year(c)

# add'l methods
year_fraction(c::Union{Actual360, Thirty360, Actual365}, d_start::Date, d_end::Date, ::Date, ::Date) = year_fraction(c, d_start, d_end)
