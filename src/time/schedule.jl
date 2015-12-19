abstract DateGenerationRule
type DateGenerationBackwards <: DateGenerationRule end
type DateGenerationForwards <: DateGenerationRule end

type Schedule
  effectiveDate::Date
  terminationDate::Date
  tenor::TenorPeriod
  convention::BusinessDayConvention
  termDateConvention::BusinessDayConvention
  rule::DateGenerationRule
  endOfMonth::Bool
  dates::Vector{Date}

  function Schedule(effectiveDate::Date, terminationDate::Date, tenor::TenorPeriod, convention::BusinessDayConvention, termDateConvention::BusinessDayConvention,
    rule::DateGenerationRule, endOfMonth::Bool)
    # dt = effectiveDate
    # num_dates = 1
    #
    # while dt < terminationDate
    #   dt += tenor.period
    #   num_dates += 1
    # end
    #
    # dates = Vector{Date}(num_dates)
    #
    # dates[1] = effectiveDate
    # dates[end] = terminationDate
    #
    # dt = effectiveDate + tenor.period
    # i = 2
    # while dt < terminationDate
    #   dates[i] = dt
    #   dt += tenor.period
    #   i += 1
    # end

    # this way is about 5-10 microseconds faster for semiannual freq over 25 years
    dates = Vector{Date}()
    dt = effectiveDate
    push!(dates, dt)
    dt += tenor.period
    while dt < terminationDate
      push!(dates, dt)
      dt += tenor.period
    end

    push!(dates, terminationDate)

    new(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates)
  end
end
