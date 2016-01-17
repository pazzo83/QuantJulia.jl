include("src/QuantJulia.jl")
using QuantJulia

function generate_flatforward_ts{C <: QuantJulia.Time.BusinessCalendar}(cal::C, settlementDate::Date)
  flat_rate = Quote(0.04875825)

  ffts = FlatForwardTermStructure(settlementDate, cal, flat_rate, QuantJulia.Time.Actual365())

  return ffts
end

function main()
  cal = QuantJulia.Time.TargetCalendar()
  settlementDate = Date(2002, 2, 19)
  todays_date = Date(2002, 2, 15)
  set_eval_date!(settings, todays_date)

  const numRows = 5
  const numCols = 5
  const swaptionVols = [0.1490, 0.1340, 0.1228, 0.1189, 0.1148, 0.1290, 0.1201, 0.1146, 0.1108, 0.1040, 0.1149, 0.1112, 0.1070, 0.1010, 0.0957, 0.1047, 0.1021, 0.0980, 0.0951, 0.1270, 0.1000, 0.0950, 0.0900, 0.1230, 0.1160]
  const swaptionMats = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]
  const swaptionLengths = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]

  # flat yield term strucutre implying 1x5 swap at 5%
  rhTermStructure = generate_flatforward_ts(cal, settlementDate)

  # Define the ATM/OTM/ITM swaps
  fixedLegFrequency = QuantJulia.Time.Annual()
  fixedLegConvention = QuantJulia.Time.Unadjusted()
  floatingLegConvention = QuantJulia.Time.ModifiedFollowing()
  fixedLegDayCounter = QuantJulia.Time.EuroThirty360()
  floatingLegFrequency = QuantJulia.Time.Semiannual()

  swapType = Payer()
  dummyFixedRate = 0.03
  indexSixMonths = euribor_index(QuantJulia.Time.TenorPeriod(Dates.Month(6)), rhTermStructure)

  startDate = QuantJulia.Time.advance(Dates.Year(1), cal, settlementDate, floatingLegConvention)
  maturity = QuantJulia.Time.advance(Dates.Year(5), cal, startDate, floatingLegConvention)

  fixedSchedule = QuantJulia.Time.Schedule(startDate, maturity, QuantJulia.Time.TenorPeriod(fixedLegFrequency), fixedLegConvention, fixedLegConvention, QuantJulia.Time.DateGenerationForwards(), false, cal)
  floatSchedule = QuantJulia.Time.Schedule(startDate, maturity, QuantJulia.Time.TenorPeriod(floatingLegFrequency), floatingLegConvention, floatingLegConvention, QuantJulia.Time.DateGenerationForwards(), false, cal)

  swap = VanillaSwap(swapType, 1000.0, fixedSchedule, dummyFixedRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

  fixedATMRate = fair_rate(swap)
  fixedOTMRate = fixedATMRate * 1.2
  fixedITMRate = fixedATMRate * 0.8

  times = zeros(0)

  for i = 1:numRows
    j = numCols - (i - 1)
    k = (i - 1) * numCols + j

    sh = SwaptionHelper(swaptionMats[i], swaptionLengths[j], Quote(swaptionVols[k]), indexSixMonths, indexSixMonths.tenor, indexSixMonths.dc, indexSixMonths.dc, rhTermStructure)

    times = add_times_to!(sh, times)
  end

  tg = QuantJulia.Time.TimeGrid(times, 30)

  # models
  modelG2 = G2(rhTermStructure)

  return modelG2
end
