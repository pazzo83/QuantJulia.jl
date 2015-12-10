# Math module
module Math

# interpolation.jl
export Interpolation, LogInterpolation, update!, locate, initialize!, value

include("interpolation.jl")

# solvers.jl
export Solver1D, BrentSolver, NewtonSolver, FiniteDifferenceNewtonSafe, solve

include("solvers.jl")

# optimization.jl
export CostFunction, Constraint, NoConstraint

include("optimization.jl")

# misc functions - prob put in own file
function divide_array_by_self!(a::Vector{Float64}, x::Number)
  for i = 1:length(a)
    a[i] = a[i] / x
  end

  return a
end

export divide_array_by_self!

end
