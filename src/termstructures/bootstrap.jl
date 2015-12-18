using QuantJulia.Math

# bootstrapping traits
const avg_rate = 0.05
const max_rate = 1.0

# Discount bootstrap trait
type Discount <: BootstrapTrait end

initial_value(::Discount) = 1.0
max_iterations(::Discount) = 100

function guess(::Discount, i::Int64, ts::TermStructure, valid::Bool)
  if valid
    # return previous iteration value
    return ts.data[i]
  end

  if i == 2
    # first pillar
    return 1.0 / (1.0 + avg_rate * ts.times[2])
  end

  r = -log(ts.data[i - 1]) / ts.times[i - 1]
  return exp(-r * ts.times[i])
end

function min_value_after(::Discount, i::Int64, ts::TermStructure, valid::Bool)
  if valid
    return ts.data[end] / 2.0
  end

  dt = ts.times[i] - ts.times[i - 1]
  return ts.data[i - 1] * exp( -max_rate * dt)
end

function max_value_after(::Discount, i::Int64, ts::TermStructure)
  return ts.data[i - 1]
end

function update_guess!(::Discount, i::Int64, ts::TermStructure, discount::Float64)
  ts.data[i] = discount
  return ts
end

# BOOTSTRAPPING
type IterativeBootstrap <: Bootstrap end

# returns initial bootstrap state to Term Structure
function initialize(::IterativeBootstrap, ts::TermStructure)
  # get the intial data based on trait
  data_initial = initial_value(ts.trait)
  n = length(ts.instruments) + 1

  # initialize data
  if data_initial == 1.0
    ts.data = ones(n)
  elseif data_initial == 0.0
    ts.data = zeros(n)
  else
    ts.data = fill(data_initial, n)
  end

  # build times and error vectors (which have the functions for the solver)
  ts.times[1] = time_from_reference(ts, ts.reference_date)
  for i = 2:n
    ts.times[i] = time_from_reference(ts, ts.instruments[i - 1].maturityDate)
    ts.errors[i] = bootstrap_error(ts, i, ts.instruments[i - 1])
  end

  # initialize interpolation
  QuantJulia.Math.initialize!(ts.interp, ts.times, ts.data)
end

function calculate!(::Bootstrap, ts::TermStructure, solver::Solver1D, first_solver::Solver1D)
  max_iter = max_iterations(ts.trait)
  valid_data = ts.validCurve

  iterations = 0
  # if we get through this loop, we haven't converged
  while iterations < max_iter
    prev_data = copy(ts.data) # need actual copy, not pointer, to check later for convergence
    for i = 2:length(ts.data)

      # bracket root and calculate guess
      min = min_value_after(ts.trait, i, ts, valid_data)
      max = max_value_after(ts.trait, i, ts)
      g = guess(ts.trait, i, ts, valid_data)

      # adjust if needed
      if g >= max
        g = max - (max - min) / 5.0
      elseif (g <= min)
        g = min + (max - min) / 5.0
      end

      if !valid_data
        update!(ts.interp, i)
      end

      # put this in a try / catch
      if !valid_data
        # use first solver
        root = solve(first_solver, ts.errors[i], ts.accuracy, g, min, max)
      else
        root = solve(solver, ts.errors[i], ts.accuracy, g, min, max)
      end
    end

    # let's check for convergence
    change = abs(ts.data[2] - prev_data[2])
    for i=3:length(ts.data)
      change = max(change, abs(ts.data[i] - prev_data[i]))
    end
    if change <= ts.accuracy
      break # bye
    end

    iterations += 1

    valid_data = true
  end
  ts.validCurve = true

  return ts
end

function bootstrap_error(ts::TermStructure, i::Int64, inst::Instrument)
  function bootstrap_error_inner(g::Float64)
    # update trait
    update_guess!(ts.trait, i, ts, g)
    update!(ts.interp, i, g)
    return quote_error(ts, inst)
  end

  return bootstrap_error_inner
end

quote_error(ts::TermStructure, inst::Instrument) = value(inst) - implied_quote(inst, ts, inst.pricing_engine)
