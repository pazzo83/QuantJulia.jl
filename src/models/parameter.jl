type ConstantParameter <: Parameter
  data::Vector{Float64}
  constraint::Constraint
end

type G2FittingParameter{T <: TermStructure} <: Parameter
  a::Float64
  sigma::Float64
  b::Float64
  eta::Float64
  rho::Float64
  ts::T
end

test_params(c::ConstantParameter, params::Vector{Float64}) = QuantJulia.Math.test(c.constraint, params)
