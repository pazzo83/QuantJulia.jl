# Curves

abstract Curve

type PiecewiseLogCubicDiscountCurve <: Curve
  calculation_date::Date
  jump_dates::Vector{JumpDate}
  dc::DayCounter
end
