# Curves

using QuantJulia.Math, QuantJulia.Time

type PiecewiseYieldCurve{I, T, B} <: InterpolatedCurve{I, T, B}
  reference_date::Date
  instruments::Vector{Instrument}
  dc::DayCount
  interp::I
  trait::T
  accuracy::Float64
  boot::B
  times::Vector{Float64}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve{I, T, B}(reference_date::Date, instruments::Vector{Instrument}, dc::DayCount, interp::I, trait::T, accuracy::Float64, boot::B)
  # get the initial length of instruments
  n = len(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve{I, T, B}(reference_date,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n),
                            Vector{Float64}(n),
                            Vector{Function}(n),
                            false)

  # initialize the bootstrapping
  initialize(boot, pyc)

  return pyc
end


max_date(curve::InterpolatedCurve) = curve.dates[end]

function discount_impl(curve::InterpolatedCurve, t::Float64)
  if t < curve.times[end]
    return value(curve.interp, t)
  end

  # do flat fwd extrapolation
end

function calculate!(curve::InterpolatedCurve)
  calculate!(curve.boot)

  return curve
end
