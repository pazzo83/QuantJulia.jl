type FittingMethodCommons{T}
  solution::Vector{Float64}
  guessSolution::Vector{Float64}
  numberOfIterations::Integer
  minimumCostValue::Float64
  weights::Vector{T}
  costFunction::FittingCost
end

function FittingMethodCommons(size::Integer, gsize::Integer)
  solution = zeros(size)
  guessSolution = zeros(gsize)
  numberOfIterations = 0
  minimumCostValue = 0.0
  weights = zeros(size)
  curve = NullCurve()
  costFunction = FittingCost(size, curve)

  return FittingMethodCommons(solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction)
end

type ExponentialSplinesFitting <: FittingMethod
  constrainAtZero::Bool
  size::Integer
  commons::FittingMethodCommons
end

function ExponentialSplinesFitting(constrainAtZero::Bool, size::Integer)
  if constrainAtZero
    gsize = 9
  else
    gsize = 10
  end

  commons = FittingMethodCommons(size, gsize)

  return ExponentialSplinesFitting(constrainAtZero, gsize, commons)
end

type SimplePolynomialFitting <: FittingMethod
  constrainAtZero::Bool
  degree::Integer
  size::Integer
  commons::FittingMethodCommons
end

function SimplePolynomialFitting(constrainAtZero::Bool, degree::Integer, size::Integer)
  if constrainAtZero
    gsize = degree
  else
    gsize = degree + 1
  end

  commons = FittingMethodCommons(size, gsize)

  return SimplePolynomialFitting(constrainAtZero, degree, gsize, commons)
end

type NelsonSiegelFitting <: FittingMethod
  constrainAtZero::Bool
  size::Integer
  commons::FittingMethodCommons
end

function NelsonSiegelFitting(size::Integer)
  constrainAtZero = true
  gsize = 4

  commons = FittingMethodCommons(size, gsize)

  return NelsonSiegelFitting(constrainAtZero, gsize, commons)
end

type SvenssonFitting <: FittingMethod
  constrainAtZero::Bool
  size::Integer
  commons::FittingMethodCommons
end

function SvenssonFitting(size::Integer)
  constrainAtZero = true
  gsize = 6

  commons = FittingMethodCommons(size, gsize)

  return SvenssonFitting(constrainAtZero, gsize, commons)
end

# function ExponentialSplinesFitting(constrainAtZero::Bool, size::Integer)
#   solution = zeros(size)
#   if constrainAtZero
#     gsize = 9
#   else
#     gsize = 10
#   end
#   guessSolution = zeros(gsize)
#   numberOfIterations = 0
#   minimumCostValue = 0.0
#   weights = zeros(size)
#   curve = NullCurve()
#   costFunction = FittingCost(size, curve)
#
#   return ExponentialSplinesFitting(constrainAtZero, gsize,
#           FittingMethodCommons(solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction))
# end

guess_size(fitting::ExponentialSplinesFitting) = fitting.constrainAtZero ? 9 : 10
guess_size(fitting::SimplePolynomialFitting) = fitting.constrainAtZero ? fitting.degree : fitting.degree + 1

# Discount functions
function discount_function{T}(method::ExponentialSplinesFitting, x::Vector{T}, t::Float64)
  d = 0.0
  N = guess_size(method)
  kappa = x[N]
  coeff = 0.0
  if !method.constrainAtZero
    for i = 1:N
      @inbounds d += x[i] * exp(-kappa * (i) * t)
    end
  else
    for i = 1:N - 1
      @inbounds d += x[i]  * exp(-kappa * (i + 1) * t)
      @inbounds coeff += x[i]
    end
    coeff = 1.0 - coeff
    d += coeff * exp(-kappa * t)
  end
  return d
end

function discount_function{T}(method::SimplePolynomialFitting, x::Vector{T}, t::Float64)
  d = 0.0
  N = method.size

  if !method.constrainAtZero
    for i = 1:N
      @inbounds d += x[i] * get_polynomial(BernsteinPolynomial(), i-1, i-1, t)
    end
  else
    d = 1.0
    for i = 1:N
      @inbounds d += x[i] * get_polynomial(BernsteinPolynomial(), i, i, t)
    end
  end

  return d
end

function discount_function{T}(method::NelsonSiegelFitting, x::Vector{T}, t::Float64)
  kappa = x[method.size]
  zero_rate = x[1] + (x[2] + x[3]) * (1.0 - exp(-kappa * t)) / ((kappa + QuantJulia.Math.EPS_VAL) * (t + QuantJulia.Math.EPS_VAL)) - (x[3]) * exp(-kappa * t)
  d = exp(-zero_rate * t)

  return d
end

function discount_function{T}(method::SvenssonFitting, x::Vector{T}, t::Float64)
  kappa = x[method.size - 1]
  kappa_1 = x[method.size]

  zero_rate = x[1] + (x[2] + x[3]) * (1.0 - exp(-kappa * t)) / ((kappa + QuantJulia.Math.EPS_VAL) * (t + QuantJulia.Math.EPS_VAL)) -
              (x[3]) * exp(-kappa * t) + x[4] * (((1.0 - exp(-kappa_1 * t)) / ((kappa_1 + QuantJulia.Math.EPS_VAL) * (t + QuantJulia.Math.EPS_VAL))) - exp(-kappa_1 * t))

  d = exp(-zero_rate * t)
  return d
end
