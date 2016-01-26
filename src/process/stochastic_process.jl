type OrnsteinUhlenbeckProcess <: StochasticProcess1D
  speed::Float64
  vol::Float64
  x0::Float64
  level::Float64
end

OrnsteinUhlenbeckProcess(speed::Float64, vol::Float64, x0::Float64 = 0.0, level::Float64 = 0.0) = OrnsteinUhlenbeckProcess(speed, vol, x0, level)

function variance(process::OrnsteinUhlenbeckProcess, ::Float64, ::Float64, dt::Float64)
  v = process.vol
  if process.speed < sqrt(eps())
    # algebraic limits for small speed
    return v * v * dt
  else
    return 0.5 * v * v / process.speed * (1.0 - exp(-2.0 * process.speed * dt))
  end
end
