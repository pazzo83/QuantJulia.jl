# Math module
module Math

# types for derivatives, etc
abstract FunctionType

type Derivative <: FunctionType end

# misc functions - prob put in own file
function divide_array_by_self!(a::Vector{Float64}, x::Number)
  for i = 1:length(a)
    a[i] = a[i] / x
  end

  return a
end

function multiply_array_by_self!(a::Vector{Float64}, x::Number)
  for i = 1:length(a)
    a[i] = a[i] * x
  end

  return a
end

export divide_array_by_self!, multiply_array_by_self!, FunctionType, Derivative

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
