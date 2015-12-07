# Curves

using QuantJulia.Math, QuantJulia.Time

type PiecewiseYieldCurve{I} <: InterpolatedCurve{I}
  reference_date::Date
  instruments::Vector{I}
  dc::DayCount
  interp::Interpolation
  trait::BootstrapTrait
  accuracy::Float64
  boot::Bootstrap
  times::Vector{Float64}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve{I}(reference_date::Date, instruments::Vector{I}, dc::DayCount, interp::Interpolation, trait::BootstrapTrait,
                                accuracy::Float64, boot::Bootstrap)
  # get the initial length of instruments
  n = length(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve(reference_date,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n + 1),
                            Vector{Float64}(n + 1),
                            Vector{Function}(n + 1),
                            false)

  # initialize the bootstrapping
  initialize(boot, pyc)

  return pyc
end


max_date(curve::InterpolatedCurve) = curve.dates[end]

function discount(curve::InterpolatedCurve, t::Float64)
  return discount_impl(curve, t)
end

function discount_impl(curve::InterpolatedCurve, t::Float64)
  if t <= curve.times[end]
    return QuantJulia.Math.value(curve.interp, t)
  end

  println("outside")
  println(t)
  println(curve.times)

  # do flat fwd extrapolation
end

function calculate!(curve::InterpolatedCurve)
  calculate!(curve.boot)

  return curve
end
