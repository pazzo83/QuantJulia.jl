# Math module
module Math

# interpolation.jl
export Interpolation, LogInterpolation, update!, locate, initialize!, value

include("interpolation.jl")

# solvers.jl
export Solver1D, BrentSolver, solve

include("solvers.jl")

end
