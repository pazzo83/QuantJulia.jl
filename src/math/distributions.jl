# Thank you Distributions.jl!

using Distributions

# normal distribution methods
function distribution_derivative(w::Normal, x::Float64)
  _sigma = std(w)
  xn = (x - mean(w)) / sigma
  return pdf(w, xn) / sigma
end
