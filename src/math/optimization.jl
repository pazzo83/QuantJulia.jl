abstract CostFunction

abstract Constraint

type NoConstraint <: Constraint end

abstract OptimizationMethod

type Simplex <: OptimizationMethod
  lambda::Float64
end

type Problem
  costFunction::CostFunction
  constraint::Constraint
  initialValue::Vector{Float64}
  currentValue::Vector{Float64}
  functionValue::Float64
  squaredNorm::Float64
  functionEvaluation::Integer
  gradientEvaluation::Integer
end

function Problem(costFunction::CostFunction, constraint::Constraint, initialValue::Vector{Float64})
  currentValue = Vector{Float64}(length(initialValue))
  functionValue = 0.0
  squaredNorm = 0.0
  functionEvaluation = 0
  gradientEvaluation = 0

  return Problem(costFunction, constraint, initialValue, currentValue, functionValue, squaredNorm, functionEvaluation, gradientEvaluation)
end

type EndCriteria
  maxIterations::Integer
  maxStationaryStateIterations::Integer
  rootEpsilon::Float64
  functionEpsilon::Float64
  gradientNormEpsilon::Float64
end

function update(constraint::Constraint, params::Vector{Float64}, direction::Vector{Float64}, beta::Float64)
  diff = beta
  new_params = params + (diff * direction)
  valid = test(constraint, new_params)
  icount = 0
  while !valid
    if (icount > 200)
      error("Can't update parameter vector")
    end

    diff *= 0.5
    icount += 1
    new_params = params + (diff * direction)
    valid = test(constraint, new_params)
  end

  params += diff * direction
  return params, diff
end

function minimize!(simplex::Simplex, p::Problem, end_criteria::EndCriteria)
  xtol = end_criteria.rootEpsilon
  max_stationary_state_iterations = end_criteria.maxStationaryStateIterations
  reset!(p)
  x = p.currentValue
  iter_num = 0

  # initialize the vertices of the simplex
  end_condition = false
  n = length(x)
  vertices = Vector{Vector{Float64}}(n + 1)
  fill!(vertices, x)
  direction = zeros(n)
  for i = 1:n
    direction[i] = 1.0
    vertices[i + 1], simplex.lambda = update(p.constraint, vertices[i + 1], simplex.lambda)
    # reset direction
    direction[i] = 0.0
  end

  # initialize function values at the vertices of the simplex
  values = zeros(n + 1)
  for i=1:n+1
    values[i] = value!(p, vertices[i])
  end

  # loop through looking for the minimum
  while !end_condition
    sum_array = zeros(n)


end


## Problem methods ##
function reset!(p::Problem)
  p.functionEvaluation = p.gradientEvaluation = 0
  p.functionValue = p.squaredNorm = 0.0

  return p
end

function value!(p::Problem, x::Vector{Float64})
  p.functionEvaluation += 1
  return value(p.costFunction, x)
end
