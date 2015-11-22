type TenorPeriod
  period::Dates.Period
  freq::Frequency

  function TenorPeriod(f::Frequency)
    freq = value(f)
    if freq < 0
      period = Dates.Day(0)
    elseif freq == 0
      period = Dates.Year(0)
    elseif freq <= 12
      period = Dates.Month(12 / freq)
    elseif freq <= 52
      period = Dates.Month(52/freq)
    else
      period = Dates.Day(1)
    end

    new(period, f)
  end
end
