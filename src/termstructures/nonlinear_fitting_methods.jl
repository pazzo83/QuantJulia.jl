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
