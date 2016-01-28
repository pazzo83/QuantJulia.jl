using QuantJulia.Math

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
  common::ShortRateModelCommon
end

function HullWhite{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01)
  _rate = forward_rate(ts, 0.0, 0.0, ContinuousCompounding(), NoFrequency())
  r0 = _rate.rate
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  phi  = HullWhiteFittingParameter(a, sigma, ts)

  return HullWhite(r0, a_const, sigma_const, phi, ts, privateConstraint, ShortRateModelCommon())
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

function notify_observers!(m::HullWhite)
  for obsv in m.common.observers
    update!(obsv)
  end

  return m
end

function tree(model::HullWhite, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)
  numericDynamics = HullWhiteDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{HullWhiteDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  @simd for i = 1:length(grid.times) - 1
    @inbounds discountBond = discount(model.ts, grid.times[i + 1])
    statePrices = get_state_prices!(numericTree, i)
    sz = get_size(numericTree, i)
    @inbounds dt = grid.dt[i]
    @inbounds dx = trinomial.dx[i]
    x = get_underlying(trinomial, i, 1)
    val = 0.0
    for j = 1:sz
      @inbounds val += statePrices[j] * exp(-x * dt)
      x += dx
    end
    val = log(val / discountBond) / dt
    @inbounds set_params!(phi, grid.times[i], val)
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

# type HullWhite{T <: TermStructure} <: ShortRateModel
#   r0::Float64
#   a::ConstantParameter
#   sigma::ConstantParameter
#   phi::HullWhiteFittingParameter
# end

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
