type ExponentialSplinesFitting <: FittingMethod
  constrainAtZero::Bool
  solution::Vector{Float64}
  guessSolution::Vector{Float64}
  numberOfIterations::Integer
  minimumCostValue::Float64
  weights::Vector{Float64}
  costFunction::FittingCost
end

function ExponentialSplinesFitting(constrainAtZero::Bool, size::Integer)
  solution = zeros(size)
  if constrainAtZero
    gsize = 9
  else
    gsize = 10
  end
  guessSolution = zeros(gsize)
  numberOfIterations = 0
  minimumCostValue = 0.0
  weights = zeros(size)
  curve = NullCurve()
  costFunction = FittingCost(size, curve)

  return ExponentialSplinesFitting(constrainAtZero, solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction)
end

guess_size(fitting::ExponentialSplinesFitting) = fitting.constrainAtZero ? 9 : 10

function discount_function(method::ExponentialSplinesFitting, x::Vector{Float64}, t::Float64)
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
