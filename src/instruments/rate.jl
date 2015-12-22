type DepositRate <: AbstractRate
  rate::Quote
  tenor::Base.Dates.Period
  fixingDays::Integer
  calendar::BusinessCalendar
  convention::BusinessDayConvention
  endOfMonth::Bool
  dc::DayCount
  iborIndex::IborIndex
  evaluationDate::Date
  referenceDate::Date
  earliestDate::Date
  maturityDate::Date
  fixingDate::Date
end

function DepositRate(rate::Quote, tenor::Base.Dates.Period, fixingDays::Integer, calendar::BusinessCalendar, convention::BusinessDayConvention, endOfMonth::Bool, dc::DayCount)
  ibor_index = IborIndex("no-fix", tenor, fixingDays, NullCurrency(), calendar, convention, endOfMonth, dc)
  evaluation_date = settings.evaluation_date
  reference_date = adjust(calendar, convention, evaluation_date)
  earliest_date = advance(Dates.Day(fixingDays), calendar, reference_date, convention)
  maturity_d = maturity_date(ibor_index, earliest_date)
  fix_date = fixing_date(ibor_index, earliest_date)
  return DepositRate(rate, tenor, fixingDays, calendar, convention, endOfMonth, dc, ibor_index, evaluation_date, reference_date, earliest_date, maturity_d, fix_date)
end
