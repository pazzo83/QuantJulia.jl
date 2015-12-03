# Curves

abstract Curve <: TermStructure

type PiecewiseYieldCurve <: Curve
  reference_date::Date
  instruments::Vector{Instrument}
  dc::DayCounter
  interp::Interpolation
  trait::BootstrapTrait
  accuracy::Float64
  boot::Bootstrap
  times::Vector{Int64}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve(reference_date::Date, instruments::Vector{Instrument}, dc::DayCounter, interp::Interpolation, trait::BootstrapTrait, accuracy::Float64, boot::Bootstrap)
  # get the initial length of instruments
  n = len(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve(reference_date,
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
