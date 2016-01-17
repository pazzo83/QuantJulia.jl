type G2{T <: TermStructure} <: ShortRateModel
  a::ConstantParameter
  sigma::ConstantParameter
  b::ConstantParameter
  eta::ConstantParameter
  rho::ConstantParameter
  phi::G2FittingParameter
  ts::T
end

function G2{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01, b::Float64 = 0.1, eta::Float64 = 0.01, rho::Float64 = -0.75)
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())
  b_const = ConstantParameter([b], PositiveConstraint())
  eta_const = ConstantParameter([eta], PositiveConstraint())
  rho_const = ConstantParameter([rho], BoundaryConstraint(-1.0, 1.0))

  phi = G2FittingParameter(a, sigma, b, eta, rho, ts)

  return G2(a_const, sigma_const, b_const, eta_const, rho_const, phi, ts)
end

# type HullWhite{T <: TermStructure} <: ShortRateModel
#   r0::Float64
#   a::ConstantParameter
#   sigma::ConstantParameter
#   phi::HullWhiteFittingParameter
# end
