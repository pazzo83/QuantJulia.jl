# bootstrapping traits
const avg_rate = 0.05
const max_rate = 1.0

abstract BootstrapTrait

# Discount bootstrap trait
type Discount <: BootstrapTrait end

initial_value(::Discount) = 1.0
max_iterations(::Discount) = 100

function guess(::Discount, i::Int64, boot::Bootstrap)
  if boot.validCurve
    # return previous iteration value
    return boot.data[i]
  end
  
  if i == 2
    # first pillar
    return 1.0 / (1.0 + avg_rate * boot.times[2])
  end
  
  r = -log(boot.data[i - 1]) / boot.times[i - 1]
  return exp(-r * boot.times[i])
end

function min_value_after(::Discount, i::Int64, boot::Bootstrap)
  if boot.validCurve
    return boot.data[end] / 2.0
  end
  
  dt = boot.times[i] - boot.times[i - 1]
  return boot.data[i - 1] * exp( -max_rate * dt)
end

function max_value_after(::Discount, i::Int64, boot::Bootstrap)
  return boot.data[i - 1]
end

function update_guess!(::Discount, i::Int64, boot::Boostrap, discount::Float64)
  boot.data[i] = discount
  return boot
end

# BOOTSTRAPPING
abstract Bootstrap

type IterativeBootstrap <: Bootstrap
  ts::TermStructure
  times::Vector{Int64}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
  
  function IterativeBootstrap(ts::TermStructure)
    # setup and initialize constructor
    n = len(ts.instruments)
    times = zeros(n)
    errors = Vector{Function}(n)
    data_initial = initial_value(ts.trait)
    
    # initialize data
    if data_initial == 1.0
      data = ones(n)
    elseif data_initial == 0.0
      data = zeros(n)
    else
      data = fill(data_initial, n)
    end
    
    new(ts, times, data, false)
  end
end

function initialize!(boot::Bootstrap)
  n = len(boot.times)
  for i = 1:n
    boot.times[i] = time_from_reference(boot.ts, boot.ts.instruments[i].maturityDate)
    boot.errors[i] = bootstrap_error(boot, ts.instruments[i], i)
  end
  
  return boot
end

function calculate!(boot::Bootstrap)
  max_iter = max_iterations(boot.ts.trait)
  
  iterations = 0
  # if we get through this loop, we haven't converged
  while iterations < max_iter
    for i = 2:len(boot.data)
      
      # bracket root and calculate guess
      min = min_value_after(boot.ts.trait, i, boot)
      max = max_value_after(boot.ts.trait, i, boot)
      guess = guess(boot.ts.trait, i, boot)
      
      # adjust if needed
      if guess >= max
        guess = max - (max - min) / 5.0
      elseif (guess <= min)
        guess = min + (max - min) / 5.0
      end
      
      if !boot.validCurve
        update!(boot.ts.interp, i)
      end
      
      # put this in a try / catch
      if !boot.validCurve
        # stuff
      end
    end
  end
end

function bootstrap_error(boot::Bootstrap, inst::Instrument, i::Int64)
  function bootstrap_error_inner(guess::Float64)
    # update trait
    update_guess!(boot.ts.trait, i, boot, guess)
    update!(boot.ts.interp, i)
    return quote_error(boot.ts, inst)
  end
  
  return bootstrap_error_inner
end

quote_error(ts::TermStructure, inst::Instrument) = value(inst) - implied_quote(ts, inst)
