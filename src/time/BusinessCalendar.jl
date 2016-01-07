# Business Calendars (adapted from Ito.jl and BusinessDays.jl)
using Base.Dates

abstract BusinessCalendar

abstract WesternCalendar <: BusinessCalendar
abstract OrthodoxCalendar <: BusinessCalendar

# target calendar
type TargetCalendar <: BusinessCalendar end

# US Calendars
abstract UnitedStatesCalendar <: WesternCalendar

type USSettlementCalendar <: UnitedStatesCalendar; end
type USNYSECalendar <: UnitedStatesCalendar; end
type USNERCCalendar <: UnitedStatesCalendar; end
type USGovernmentBondCalendar <: UnitedStatesCalendar; end

abstract BusinessDayConvention
type Unadjusted <: BusinessDayConvention end
type ModifiedFollowing <: BusinessDayConvention end
type Following <: BusinessDayConvention end

# easter functions
function easter_rata{I <: Integer}(y::I)

  c::Int64
	e::Int64
	p::Int64

   # Algo R only works after 1582
   if y < 1582
        # Are you using this? Send me a postcard!
        error("Year cannot be less than 1582. Provided: $(y).")
   end

	# Century
   c = div( y , 100) + 1

   # Shifted Epact
   e = mod(14 + 11*(mod(y, 19)) - div(3*c, 4) + div(5+8*c, 25), 30)

   # Adjust Epact
   if (e == 0) || ((e == 1) && ( 10 < mod(y, 19) ))
   	e += 1
   end

   # Paschal Moon
   p = Date(y, 4, 19).instant.periods.value - e

   # Easter: locate the Sunday after the Paschal Moon
   return p + 7 - mod(p, 7)
end

# Returns Date
function easter_date{I <: Integer}(y::I)
	# Compute the gregorian date for Rata Die number
     return Date(Dates.rata2datetime( easter_rata(y) ))
end

# calendar functions
function advance{C <: BusinessCalendar, B <: BusinessDayConvention}(days::Day, cal::C, dt::Date, biz_conv::B)
  n = Int(days)
  if n > 0
    while n > 0
      dt += Day(1)
      while !is_business_day(cal, dt)
        dt += Day(1)
      end
      n -= 1
    end
  else
    while (n < 0)
      dt -= Day(1)
      while !is_business_day(cal, dt)
        dt -= Day(1)
      end
      n += 1
    end
  end

  return dt
end

function advance{C <: BusinessCalendar, B <: BusinessDayConvention}(time_period::Union{Week, Month, Year}, cal::C, dt::Date, biz_conv::B)
  dt += time_period
  return adjust(cal, biz_conv, dt)
end


function is_business_day{C <: BusinessCalendar}(cal::C, dt::Date)
  if dayofweek(dt) in [6, 7] || is_holiday(cal, dt)
    return false
  else
    return true
  end
end

# In the United States, if a holiday falls on Saturday, it's observed on the preceding Friday.
# If it falls on Sunday, it's observed on the next Monday.
function adjustweekendholidayUS(dt::Date)
	if dayofweek(dt) == 6
		return dt - Dates.Day(1)
	end

	if dayofweek(dt) == 7
		return dt + Dates.Day(1)
	end

	return dt
end

function is_holiday(::USSettlementCalendar , dt::Date)

	const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)

	if (
			# New Year's Day
			adjustweekendholidayUS(Date(yy, 1, 1)) == dt
			||
			# New Year's Day on the previous year when 1st Jan is Saturday
			(mm == 12 &&  dd == 31 && dayofweek(dt) == Friday)
			||
			# Birthday of Martin Luther King, Jr.
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 1)
			||
			# Washington's Birthday
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 2)
			||
			# Memorial Day is the last Monday in May
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
			||
			# Independence Day
			adjustweekendholidayUS(Date(yy, 7, 4)) == dt
			||
			# Labor Day is the first Monday in September
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 9)
			||
			# Columbus Day is the second Monday in October
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 2 && mm == 10)
			||
			# Veterans Day
			adjustweekendholidayUS(Date(yy, 11, 11)) == dt
			||
			# Thanksgiving Day is the fourth Thursday in November
			(dayofweek(dt) == 4 && dayofweekofmonth(dt) == 4 && mm == 11)
			||
			# Christmas
			adjustweekendholidayUS(Date(yy, 12, 25)) == dt
		)
		return true
	end

	return false
end

function is_holiday(::USGovernmentBondCalendar, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)
	if (
			# New Year's Day
			adjustweekendholidayUS(Date(yy, 1, 1)) == dt
			||
			# New Year's Day on the previous year when 1st Jan is Saturday
			(mm == 12 &&  dd == 31 && dayofweek(dt) == Friday)
			||
			# Birthday of Martin Luther King, Jr.
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 1)
			||
			# Washington's Birthday
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 2)
			||
      # Good Friday
      easter_date(yy) - Day(2) == dt
      ||
			# Memorial Day is the last Monday in May
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
			||
			# Independence Day
			adjustweekendholidayUS(Date(yy, 7, 4)) == dt
			||
			# Labor Day is the first Monday in September
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 9)
			||
			# Columbus Day is the second Monday in October
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 2 && mm == 10)
			||
			# Veterans Day
			adjustweekendholidayUS(Date(yy, 11, 11)) == dt
			||
			# Thanksgiving Day is the fourth Thursday in November
			(dayofweek(dt) == 4 && dayofweekofmonth(dt) == 4 && mm == 11)
			||
			# Christmas
			adjustweekendholidayUS(Date(yy, 12, 25)) == dt
		)
		return true
	end

	return false
end

function is_holiday(::TargetCalendar, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)
  const easter_sun = easter_date(yy)

  if (
    # New Years Day
    (mm == 1 && dd == 1)
    ||
    # Good Friday
    easter_sun - Day(2) == dt
    ||
    # Easter Monday
    easter_sun + Day(1) == dt
    ||
    # Int'l Labor Day
    (mm == 5 && dd == 1)
    ||
    # Christmas
    (mm == 12 && dd == 25)
    ||
    # Day of Goodwill
    (mm == 12 && dd == 26)
    )
    return true
  else
    return false
  end
end

# adjustments
adjust{B <: BusinessCalendar}(::B, ::Unadjusted, d::Date) = d

function adjust{B <: BusinessCalendar}(cal::B, ::Union{ModifiedFollowing, Following}, d::Date)
  while !is_business_day(cal, d)
    d += Day(1)
  end

  return d
end
