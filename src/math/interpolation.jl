# interpolation
abstract Interpolation

# Log Linear interpolation
type LogInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  interpolator::LinearInterpolation
end

function LogInterpolation(x_vals::Vector{Float64}, y_vals::Vector{Float64})
  # build log of y values, defaulting to 0 for initial state
  n = length(x_vals)
  log_y_vals = zeros(n)
  s = zeros(n) # first derivative

  # initialize the linear interpolator
  interpolator = LinearInterpolation(x_vals, log_y_vals, s)

  return LogInterpolation(x_vals, y_vals, interpolator)
end

type LinearInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  s::Vector{Float64}
end

# Log Interpolation update
function update!(interp::LogInterpolation, idx::Int64)
  # first get the log of the y values
  for i = 1:idx
    interp.interpolator.y_vals[i] = log(interp.y_vals[i])
  end

  # use these log values to update the linear interpolator
  update!(interp.interpolator, idx)

  return interp
end

# Linear Interpolation update
function update!(interp::LinearInterpolation, idx::Int64)
  for i = 2:idx
    dx = interp.x_vals[i] - interp.x_vals[i - 1]
    s[i - 1] = (interp.y_vals[i] - interp.y_vals[i - 1]) / dx
  end

  return interp
end

# locate x
function locate(interp::Interpolation, val::Float64)
  if x < interp.x_vals[1]
    return 0
  elseif x > interp.x_vals[end - 1]
    return interp.x_vals[end] - interp.x_vals[1] - 2
  else
    return findfirst(x[1:end - 1] .> 1) - x[1] - 1
  end
end

value(interp::LogInterpolation, val::Float64) = exp(value(interp.interpolator, val))

function value(interp::LinearInterpolation, val::Float64)
  i = locate(interp, val)
  return interp.y_vals[i] + (val - interp.x_vals[i]) * interp.s[i]
end
