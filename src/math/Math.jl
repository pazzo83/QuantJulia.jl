# Math module
module Math

# types for derivatives, etc
abstract FunctionType

type Derivative <: FunctionType end

type BernsteinPolynomial end

# misc functions - prob put in own file
function divide_array_by_self!{T}(a::Vector{T}, x::Number)
  for i = 1:length(a)
    a[i] = a[i] / x
  end

  return a
end

function multiply_array_by_self!{T}(a::Vector{T}, x::Number)
  for i = 1:length(a)
    a[i] = a[i] * x
  end

  return a
end

function get_factorial(i::Integer)
  if i > 20
    return Float64(factorial(BigInt(i)))
  else
    return factorial(i)
  end
end

function get_polynomial(::BernsteinPolynomial, i::Integer, n::Integer, x::Float64)
  coeff = get_factorial(n) / (get_factorial(n-1) * get_factorial(i))

  return coeff * (x ^ i) * (1.0 - x)^(n - i)
end

export divide_array_by_self!, multiply_array_by_self!, get_factorial, get_polynomial, FunctionType, Derivative, BernsteinPolynomial

# Constants
const EPS_VAL = eps()

# interpolation.jl
export Interpolation, LogInterpolation, update!, locate, initialize!, value

include("interpolation.jl")

# solvers.jl
export Solver1D, BrentSolver, NewtonSolver, FiniteDifferenceNewtonSafe, solve

include("solvers.jl")

# optimization.jl
export CostFunction, Constraint, NoConstraint, OptimizationMethod, Simplex, Problem, EndCriteria, minimize!, minimize_2!

include("optimization.jl")

end
