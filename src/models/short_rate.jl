using QuantJulia.Math
using Distributions

type PrivateConstraint{P <: Parameter} <: Constraint
  arguments::Vector{P}
end

function QuantJulia.Math.test(c::PrivateConstraint, x::Vector{Float64})
  k = 1
  for i = 1:length(c.arguments)
    sz = length(c.arguments[i].data)
    testParams = zeros(sz)
    for j = 1:sz
      testParams[j] = x[k]
      k += 1
    end
    if !test_params(c.arguments[i], testParams)
      return false
    end
  end

  return true
end

type CalibrationFunction{M <: ShortRateModel, C <: CalibrationHelper} <: CostFunction
  model::M
  helpers::Vector{C}
  weights::Vector{Float64}
  projection::Projection
end

type SolvingFunction
  lambda::Vector{Float64}
  Bb::Vector{Float64}
end

function operator(solvFunc::SolvingFunction)
  function _inner(y::Float64)
    value_ = 1.0
    for i = 1:length(solvFunc.lambda)
      value_ -= solvFunc.lambda[i] * exp(-solvFunc.Bb[i] * y)
    end
    return value_
  end

  return _inner
end

function func_values(calibF::CalibrationFunction, params::Vector{Float64})
  set_params!(calibF.model, params)
  values = zeros(length(calibF.helpers))
  for i = 1:length(values)
    values[i] = calibration_error(calibF.helpers[i].calibCommon.calibrationErrorType, calibF.helpers[i]) * sqrt(calibF.weights[i])
  end

  return values
end

function value(calibF::CalibrationFunction, params::Vector{Float64})
  set_params!(calibF.model, params)
  _value = 0.0
  for i = 1:length(calibF.helpers)
    diff = calibration_error(calibF.helpers[i].calibCommon.calibrationErrorType, calibF.helpers[i])
    _value += diff * diff * calibF.weights[i]
  end

  return sqrt(_value)
end

# accessor methods ##
get_a{M <: ShortRateModel}(m::M) = m.a.data[1]
get_sigma{M <: ShortRateModel}(m::M) = m.sigma.data[1]
get_b{M <: ShortRateModel}(m::M) = m.b.data[1]

type G2{T <: TermStructure} <: ShortRateModel
  a::ConstantParameter
  sigma::ConstantParameter
  b::ConstantParameter
  eta::ConstantParameter
  rho::ConstantParameter
  phi::G2FittingParameter
  ts::T
  privateConstraint::PrivateConstraint
end

function G2{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01, b::Float64 = 0.1, eta::Float64 = 0.01, rho::Float64 = -0.75)
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())
  b_const = ConstantParameter([b], PositiveConstraint())
  eta_const = ConstantParameter([eta], PositiveConstraint())
  rho_const = ConstantParameter([rho], BoundaryConstraint(-1.0, 1.0))

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const, b_const, eta_const, rho_const])

  phi = G2FittingParameter(a, sigma, b, eta, rho, ts)

  return G2(a_const, sigma_const, b_const, eta_const, rho_const, phi, ts, privateConstraint)
end

get_eta(m::G2) = m.eta.data[1]
get_rho(m::G2) = m.rho.data[1]

get_params(m::G2) = Float64[get_a(m), get_sigma(m), get_b(m), get_eta(m), get_rho(m)]

generate_arguments!(m::G2) = m.phi = G2FittingParameter(get_a(m), get_sigma(m), get_b(m), get_eta(m), get_rho(m), m.ts)

function V(m::G2, t::Float64)
  expat = exp(-get_a(m) * t)
  expbt = exp(-get_b(m) * t)
  cx = get_sigma(m) / get_a(m)
  cy = get_eta(m) / get_b(m)
  valuex = cx * cx * (t + (2.0 * expat - 0.5 * expat * expat - 1.5) / get_a(m))
  valuey = cy * cy * (t + (2.0 * expbt - 0.5 * expbt * expbt - 1.5) / get_b(m))
  value = 2.0 * get_rho(m) * cx * cy * (t + (expat - 1.0) / get_a(m) + (expbt - 1.0) / get_b(m) - (expat * expbt - 1.0) / (get_a(m) + get_b(m)))

  return valuex + valuey + value
end

A(m::G2, t::Float64, T::Float64) = discount(m.ts, T) / discount(m.ts, t) * exp(0.5 * (V(m, T - t) - V(m, T) + V(m, t)))
B(m::G2, x::Float64, t::Float64) = (1.0 - exp(-x * t)) / x

type G2SwaptionPricingFunction
  a::Float64
  sigma::Float64
  b::Float64
  eta::Float64
  rho::Float64
  w::Float64
  startTime::Float64
  fixedRate::Float64
  sigmax::Float64
  sigmay::Float64
  rhoxy::Float64
  mux::Float64
  muy::Float64
  payTimes::Vector{Float64}
  A::Vector{Float64}
  Ba::Vector{Float64}
  Bb::Vector{Float64}
end

function G2SwaptionPricingFunction(model::G2, w::Float64, startTime::Float64, payTimes::Vector{Float64}, fixedRate::Float64)
  a = get_a(model)
  sigma = get_sigma(model)
  b = get_b(model)
  eta = get_eta(model)
  rho = get_rho(model)

  sigmax = sigma * sqrt(0.5 * (1.0 - exp(-2.0 * a * startTime)) / a)
  sigmay = eta * sqrt(0.5 * (1.0 - exp(-2.0 * b * startTime)) / b)
  rhoxy = rho * eta * sigma * (1.0 - exp(-(a + b) * startTime)) / ((a + b) * sigmax * sigmay)

  temp = sigma * sigma / (a * a)
  mux = -((temp + rho * sigma * eta / (a * b)) * (1.0 - exp(-a * startTime)) - 0.5 * temp * (1.0 - exp(-2.0 * a * startTime)) -
        rho * sigma * eta / (b * (a + b)) * (1.0 - exp(-(b + a) * startTime)))

  temp = eta * eta / (b * b)
  muy = -((temp + rho * sigma * eta / (a * b)) * (1.0 - exp(-b * startTime)) - 0.5 * temp * (1.0 - exp(-2.0 * b * startTime)) -
        rho * sigma * eta / (a * (a + b)) * (1.0 - exp(-(b + a) * startTime)))

  n = length(payTimes)
  A_ = zeros(n)
  Ba = zeros(n)
  Bb = zeros(n)
  for i = 1:n
    A_[i] = A(model, startTime, payTimes[i])
    Ba[i] = B(model, a, payTimes[i] - startTime)
    Bb[i] = B(model, b, payTimes[i] - startTime)
  end

  return G2SwaptionPricingFunction(a, sigma, b, eta, rho, w, startTime, fixedRate, sigmax, sigmay, rhoxy, mux, muy, payTimes, A_, Ba, Bb)
end

function operator(pricingFunc::G2SwaptionPricingFunction)
  phi = Normal()
  n = length(pricingFunc.payTimes)
  function _inner(x::Float64)
    temp = (x - pricingFunc.mux) / pricingFunc.sigmax
    txy = sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
    lambda = zeros(n)

    for i = 1:n
      tau = i == 1 ? pricingFunc.payTimes[1] - pricingFunc.startTime : pricingFunc.payTimes[i] - pricingFunc.payTimes[i - 1]
      c = i == n ? (1.0 + pricingFunc.fixedRate * tau) : pricingFunc.fixedRate * tau
      lambda[i] = c * pricingFunc.A[i] * exp(-pricingFunc.Ba[i] * x)
    end

    func = SolvingFunction(lambda, pricingFunc.Bb)
    s1d = BrentSolver(1000)

    yb = solve(s1d, operator(func), 1e-6, 0.00, -100.0, 100.0)

    h1 = (yb - pricingFunc.muy) / (pricingFunc.sigmay * txy) - pricingFunc.rhoxy * (x - pricingFunc.mux) / (pricingFunc.sigmax * txy)
    val = cdf(phi, -pricingFunc.w * h1)

    for i = 1:n
      h2 = h1 + pricingFunc.Bb[i] * pricingFunc.sigmay * sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
      kappa = -pricingFunc.Bb[i] * (pricingFunc.muy - 0.5 * txy * txy * pricingFunc.sigmay * pricingFunc.sigmay * pricingFunc.Bb[i] +
              pricingFunc.rhoxy * pricingFunc.sigmay * (x - pricingFunc.mux) / pricingFunc.sigmax)

      val -= lambda[i] * exp(kappa) * cdf(phi, -pricingFunc.w * h2)
    end

    blah = exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
    return exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
  end

  return _inner
end

## ONE FACTOR MODELS ##
type OneFactorShortRateTree{S <: ShortRateDynamics} <: ShortRateTree
  tree::TrinomialTree
  dynamics::S
  tg::TimeGrid
  treeLattice::TreeLattice1D

  function OneFactorShortRateTree{S}(tree::TrinomialTree, dynamics::S, tg::TimeGrid)
    oneFactorTree = new(tree, dynamics, tg)
    oneFactorTree.treeLattice = TreeLattice1D(tg, get_size(tree, 2), oneFactorTree)

    return oneFactorTree
  end
end

get_size{I <: Integer}(tr::OneFactorShortRateTree, i::I) = get_size(tr.tree, i)

get_state_prices!(tree::OneFactorShortRateTree, i::Int) = get_state_prices!(tree.treeLattice, i)

function discount(tr::OneFactorShortRateTree, i::Int, idx::Int)
  x = get_underlying(tr.tree, i, idx)
  r = short_rate(tr.dynamics, tr.tg.times[i], x)
  return exp(-r * tr.tg.dt[i])
end

descendant(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = descendant(tr.tree, i, idx, branch)
probability(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = probability(tr.tree, i, idx, branch)

type HullWhite{T <: TermStructure} <: ShortRateModel
  r0::Float64
  a::ConstantParameter
  sigma::ConstantParameter
  phi::HullWhiteFittingParameter
  ts::T
  privateConstraint::PrivateConstraint
end

function HullWhite{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01)
  _rate = forward_rate(ts, 0.0, 0.0, ContinuousCompounding(), NoFrequency())
  r0 = _rate.rate
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  phi  = HullWhiteFittingParameter(a, sigma, ts)

  return HullWhite(r0, a_const, sigma_const, phi, ts, privateConstraint)
end

type HullWhiteDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

HullWhiteDynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64) = HullWhiteDynamics(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::HullWhiteDynamics, t::Float64, x::Float64) = x + operator(dynamic.fitting, t)

get_params(m::HullWhite) = Float64[get_a(m), get_sigma(m)]

generate_arguments!(m::HullWhite) = m.phi = HullWhiteFittingParameter(get_a(m), get_sigma(m), m.ts)

function tree(model::HullWhite, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)
  numericDynamics = HullWhiteDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{HullWhiteDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  for i = 1:length(grid.times) - 1
    discountBond = discount(model.ts, grid.times[i + 1])
    statePrices = get_state_prices!(numericTree, i)
    sz = get_size(numericTree, i)
    dt = grid.dt[i]
    dx = trinomial.dx[i]
    x = get_underlying(trinomial, i, 1)
    val = 0.0
    for j = 1:sz
      val += statePrices[j] * exp(-x * dt)
      x += dx
    end
    val = log(val / discountBond) / dt
    set_params!(phi, grid.times[i], val)
  end

  return numericTree
end

type RStarFinder{M <: ShortRateModel}
  model::M
  strike::Float64
  maturity::Float64
  valueTime::Float64
  fixedPayTimes::Vector{Float64}
  amounts::Vector{Float64}
end

function operator(rsf::RStarFinder)
  function _inner(x::Float64)
    _value = rsf.strike
    _B = discount_bond(rsf.model, rsf.maturity, rsf.valueTime, x)
    sz = length(rsf.fixedPayTimes)
    for i = 1:sz
      dbVal = discount_bond(rsf.model, rsf.maturity, rsf.fixedPayTimes[i], x) / _B
      _value -= rsf.amounts[i] * dbVal
    end

    return _value
  end

  return _inner
end

function B(model::HullWhite, t::Float64, T::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    return T - t
  else
    return (1.0 - exp(-_a * (T - t))) / _a
  end
end

function A(model::HullWhite, t::Float64, T::Float64)
  discount1 = discount(model.ts, t)
  discount2 = discount(model.ts, T)

  forward = forward_rate(model.ts, t, t, ContinuousCompounding(), NoFrequency())
  temp = get_sigma(model) * B(model, t, T)
  val = B(model, t, T) * forward.rate - 0.25 * temp * temp * B(model, 0.0, 2.0* t)

  return exp(val) * discount2 / discount1
end

function set_params!{M <: ShortRateModel}(model::M, params::Vector{Float64})
  paramCount = 1
  args = model.privateConstraint.arguments
  for i = 1:length(args)
    arg = model.privateConstraint.arguments[i]
    for j = 1:length(arg.data)
       arg.data[j] = params[paramCount]
       paramCount += 1
     end
   end
   generate_arguments!(model)
   # model.phi = G2FittingParameter(get_a(model), get_sigma(model), get_b(model), get_eta(model), get_rho(model), model.ts)
   return model
 end

# type HullWhite{T <: TermStructure} <: ShortRateModel
#   r0::Float64
#   a::ConstantParameter
#   sigma::ConstantParameter
#   phi::HullWhiteFittingParameter
# end

function calibrate!{M <: ShortRateModel, C <: CalibrationHelper, O <: OptimizationMethod}(model::M, instruments::Vector{C}, method::O, endCriteria::EndCriteria,
                    constraint::Constraint = model.privateConstraint, weights::Vector{Float64} = ones(length(instruments)), fixParams::BitArray = BitArray(0))
  prms = get_params(model)
  all = falses(length(prms))
  proj = Projection(prms, length(fixParams) > 0 ? fixParams : all)
  calibFunc = CalibrationFunction(model, instruments, weights, proj)
  pc = ProjectedConstraint(constraint, proj)
  prob = Problem(calibFunc, pc, project(proj, prms))

  # minimization
  minimize!(method, prob, endCriteria)
  res = prob.currentValue
  set_params!(model, include_params(proj, res))

  return model
end

function gen_swaption{I <: Integer}(model::G2, swaption::Swaption, fixedRate::Float64, range::Float64, intervals::I)
  settlement = reference_date(model.ts)
  dc = model.ts.dc
  startTime = year_fraction(dc, settlement, swaption.swap.args.floatingResetDates[1])
  w = swaption.swap.payer[2]

  fixedPayTimes = zeros(length(swaption.swap.args.fixedPayDates))
  for i = 1:length(fixedPayTimes)
    fixedPayTimes[i] = year_fraction(dc, settlement, swaption.swap.args.fixedPayDates[i])
  end
  func = G2SwaptionPricingFunction(model, w, startTime, fixedPayTimes, fixedRate)

  upper = func.mux + range * func.sigmax
  lower = func.mux - range * func.sigmax

  integrator = SegmentIntegral(intervals)
  return swaption.swap.nominal * w * discount(model.ts, startTime) * QuantJulia.Math.operator(integrator, QuantJulia.operator(func), lower, upper)
end

discount_bond{M <: ShortRateModel}(model::M, tNow::Float64, maturity::Float64, _rate::Float64) = A(model, tNow, maturity) * exp(-B(model, tNow, maturity) * _rate)

function discount_bond_option{O <: OptionType}(model::HullWhite, optionType::O, strike::Float64, maturity::Float64, bondStart::Float64, bondMaturity::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    v = get_sigma(model) * B(model, bondStart, bondMaturity) * sqrt(maturity)
  else
    v = get_sigma(model) / (_a * sqrt(2.0 * _a)) * sqrt(exp(-2.0 * _a * (bondStart - maturity)) - exp(-2.0 * _a * bondStart) - 2.0 * (exp(-_a * (bondStart + bondMaturity - 2.0 * maturity))
        - exp(-_a * (bondStart + bondMaturity))) + exp(-2.0 * _a * (bondMaturity - maturity)) - exp(-2.0 * _a * bondMaturity))
  end

   f = discount(model.ts, bondMaturity)
   k = discount(model.ts, bondStart) * strike

   return black_formula(optionType, k, f, v)
 end
