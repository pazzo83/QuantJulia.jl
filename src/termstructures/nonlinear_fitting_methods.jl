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
  solution = Vector{Float64}(size)
  guessSolution = Vector{Float64}()
  numberOfIterations = 0
  minimumCostValue = 0.0
  weights = Vector{Float64}(size)
  curve = NullCurve()
  costFunction = FittingCost(size, curve)

  return ExponentialSplinesFitting(constrainAtZero, solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction, curve)
end

guess_size(fitting::ExponentialSplinesFitting) = fitting.constrainAtZero ? 9 : 10

function discount_function(method::ExponentialSplinesFitting, x::Vector{Float64}, t::Float64)
  d = 0.0
end
