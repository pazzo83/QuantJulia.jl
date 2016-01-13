type AmericanExercise <: EarlyExercise
  earliestDate::Date
  latestDate::Date
end

type BermudanExercise <: EarlyExercise
  dates::Vector{Date}
end

type EuropeanExercise <: Exercise
  date::Date
end
