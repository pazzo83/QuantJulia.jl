using QuantJulia.Math

# bootstrapping traits
const avg_rate = 0.05
const max_rate = 1.0

# Discount bootstrap trait
type Discount <: BootstrapTrait end

initial_value(::Discount) = 1.0
max_iterations(::Discount) = 100

function guess{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, valid::Bool)
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

function min_value_after{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, valid::Bool)
  if valid
    return ts.data[end] / 2.0
  end

  dt = ts.times[i] - ts.times[i - 1]
  return ts.data[i - 1] * exp( -max_rate * dt)
end

function max_value_after{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T)
  return ts.data[i - 1]
end

function update_guess!{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, discount::Float64)
  ts.data[i] = discount
  return ts
end

get_pricing_engine{Y}(::Discount, yts::Y) = DiscountingBondEngine(yts)

apply_termstructure!{R <: AbstractRate, T <: TermStructure}(rate::R, ts::T) = rate.iborIndex.ts = ts
apply_termstructure!{B <: Bond, T <: TermStructure}(b::B, ts::T) = b.pricingEngine.yts = ts
function apply_termstructure!{S <: Swap, T <: TermStructure}(s::S, ts::T)
  s.iborIndex.ts = ts
  s.pricingEngine.yts = ts

  return s
end

# BOOTSTRAPPING
type IterativeBootstrap <: Bootstrap end

# returns initial bootstrap state to Term Structure
function initialize{T <: TermStructure}(::IterativeBootstrap, ts::T)
  # get the intial data based on trait
  data_initial = initial_value(ts.trait)
  n = length(ts.instruments) + 1

  # get a new pricing engine based on type
  # pe = get_pricing_engine(ts.trait, ts)

  # initialize data
  if data_initial == 1.0
    ts.data = ones(n)
  elseif data_initial == 0.0
    ts.data = zeros(n)
  else
    ts.data = fill(data_initial, n)
  end

  # build times and error vectors (which have the functions for the solver)
  ts.times[1] = time_from_reference(ts, ts.referenceDate)
  for i = 2:n
    @inbounds ts.times[i] = time_from_reference(ts, maturity_date(ts.instruments[i - 1]))
    @inbounds ts.errors[i] = bootstrap_error(i, ts.instruments[i - 1], ts)
    # set yield term Structure
    @inbounds apply_termstructure!(ts.instruments[i - 1], ts)
  end

  # initialize interpolation
  QuantJulia.Math.initialize!(ts.interp, ts.times, ts.data)
end

function calculate!{T <: TermStructure, S <: Solver1D, SS <: Solver1D}(::Bootstrap, ts::T, solver::S, first_solver::SS)
  max_iter = max_iterations(ts.trait)
  valid_data = ts.validCurve
  ts.calculated = true

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

function bootstrap_error{I <: Integer, T <: Instrument, Y <: YieldTermStructure}(i::I, inst::T, ts::Y)
  function bootstrap_error_inner(g::Float64)
    # update trait
    update_guess!(ts.trait, i, ts, g)
    update!(ts.interp, i, g)
    return quote_error(inst)
  end

  return bootstrap_error_inner
end

quote_error{B <: Bond}(inst::B) = QuantJulia.value(inst) - implied_quote(inst, true) # recalculate
quote_error{R <: AbstractRate}(rate::R) = rate.rate.value - implied_quote(rate, rate.iborIndex.ts)
quote_error{S <: Swap}(swap::S) = swap.rate.value - implied_quote(swap)
