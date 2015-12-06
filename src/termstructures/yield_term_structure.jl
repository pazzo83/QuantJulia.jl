# main term structures
using QuantJulia.Time

# core TermStructure methods
function check_range(ts::TermStructure, date::Date)
  date < ts.reference_date && "Date $date before reference_date $(ts.reference_date)"
end

function check_range(ts::TermStructure, time_frac::Float64)
  time_frac < 0.0 && "Negative time given: $time"
end

max_date(ts) = ts.reference_date + Date.Year(100)

time_from_reference(ts::TermStructure, date::Date) = year_fraction(ts.dc, ts.reference_date, date)

type JumpDate
  ts_quote::Quote
  ts_date::Date
end

type JumpTime
  ts_quote::Quote
  ts_time::Float64
end

discount(yts::YieldTermStructure, date::Date) = discount(yts, time_from_reference(yts, date))

function discount(yts::YieldTermStructure, time_frac::Float64)
  disc = discount_impl(yts, time_frac)
  if length(yts.jump_times) == 0
    return disc
  end

  jump_effect = 1.0
  for jump in yts.jump_times
    if jump.ts_time > 0.0 && jump.ts_time < time
      if jump.ts_quote.value > 0.0 && jump.ts_quote.value <= 1.0
        jump_effect *= jump.ts_quote.value
      end
    end
  end

  return jump_effect * disc
end

function zero_rate(yts::YieldTermStructure, date::Date, dc::DayCount, comp::CompoundingType, freq::Frequency)
  if date == yts.reference_date
    return implied_rate(1.0 / discount(yts, 0.0001), dc, comp, 0.0001, freq)
  else
    return implied_rate(1.0 / discount(yts, date), dc, comp, date, freq)
  end
end

function zero_rate(yts::YieldTermStructure, time_frac::Float64, comp::CompoundingType, freq::Frequency)
  t = time_frac == 0.0 ? 0.0001 : time_frac
  return implied_rate(1.0 / discount(yts, t), comp, t, freq)
end

function forward_rate(yts::YieldTermStructure, date1::Date, date2::Date, dc::DayCount, comp::CompoundingType, freq::Frequency)
  if date1 == date2
    t1 = max(time_from_reference(yts, date1) - 0.0001 / 2.0, 0.0)
    t2 = t1 + 0.0001
    return implied_rate(discount(yts, t1) / discount(yts, d2), dc, comp, 0.0001, freq)
  elseif date1 < date2
    return implied_rate(discount(yts, date1) / discount(yts, date2), dc, comp, date1, date2, freq)
  else
    error("Forward start date must be before forward end date")
  end
end

forward_rate(yts::YieldTermStructure, date::Date, period::Integer, dc::DayCount, comp::CompoundingType, freq::Frequency) = forward_rate(yts, date, date + Dates.Day(period), dc, comp, freq)

function forward_rate(yts::YieldTermStructure, time1::Float64, time2::Float64, comp::CompoundingType, freq::Frequency)
  if time1 == time2
    t1 = max(time1 - 0.0001 / 2.0, 0.0)
    t2 = t1 + 0.0001
    interval, compound = (t2 - t1, discount(yts, t1) / discount(yts, t2))
  else
    interval, compound = (time2 - time1, discount(yts, time1) / discount(yts, time2))
  end

  return implied_rate(compound, dc, comp, interval, freq)
end

## FlatForwardTermStructure
type FlatForwardTermStructure <: YieldTermStructure
  settlement_days::Integer
  reference_date::Date
  calendar::BusinessCalendar
  forward::Quote
  dc::DayCount
  comp::CompoundingType
  freq::Frequency
  rate::InterestRate

  function FlatForwardTermStructure(settlement_days::Integer, reference_date::Date, calendar::BusinessCalendar, forward::Quote, dc::DayCount, comp::CompoundingType, freq::Frequency)
    rate = InterestRate(forward.value, dc, comp, freq)
    new(settlement_days, reference_date, calendar, forward, dc, comp, freq, rate)
  end
end

discount(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)

# discount_impl(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)
