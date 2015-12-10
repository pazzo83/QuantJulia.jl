# Curves

using QuantJulia.Math, QuantJulia.Time

type NullCurve <: Curve end

type PiecewiseYieldCurve{I} <: InterpolatedCurve{I}
  reference_date::Date
  instruments::Vector{I}
  dc::DayCount
  interp::Interpolation
  trait::BootstrapTrait
  accuracy::Float64
  boot::Bootstrap
  times::Vector{Float64}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve{I}(reference_date::Date, instruments::Vector{I}, dc::DayCount, interp::Interpolation, trait::BootstrapTrait,
                                accuracy::Float64, boot::Bootstrap)
  # get the initial length of instruments
  n = length(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve(reference_date,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n + 1),
                            Vector{Float64}(n + 1),
                            Vector{Function}(n + 1),
                            false)

  # initialize the bootstrapping
  initialize(boot, pyc)

  return pyc
end

type FittedBondDiscountCurve{B} <: Curve
  settlement_days::Integer
  calendar::BusinessCalendar
  bonds::Vector{B}
  dc::DayCount
  fittingMethod::FittingMethod
  accuracy::Float64
  maxEvaluations::Integer
  guessArray::Vector{Float64}
  simplexLambda::Float64
end

function FittedBondDiscountCurve{B}(settlement_days::Integer, calendar::Calendar, bonds::Vector{B}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64=1e-10,
                                    maxEvaluations=10000, simplexLambda::Float64=1.0)
  n = length(bonds)
  guessArray = Vector{Float64}(n)

  fbdc = FittedBondDiscountCurve(settlement_days, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, guessArray, simplexLambda)
  fbdc.fittingMethod.costFunction.curve = fbdc

  return fbdc
end

type FittingCost <: CostFunction
  value::Vector{Float64}
  values::Vector{Float64}
  firstCashFlow::Vector{Float64}
  curve::Curve
end

function FittingCost(size::Integer, curve::Curve)
  value = Vector{Float64}()
  values = Vector{Float64}()
  firstCashFlow = zeros(size)

  return FittingCost(value, values, firstCashFlow, curve)
end

# Interpolated Curve methods #
max_date(curve::InterpolatedCurve) = curve.dates[end]

function discount(curve::InterpolatedCurve, t::Float64)
  return discount_impl(curve, t)
end

function discount_impl(curve::InterpolatedCurve, t::Float64)
  if t <= curve.times[end]
    return QuantJulia.Math.value(curve.interp, t)
  end

  # println("outside")
  # println(t)
  # println(curve.times)

  # do flat fwd extrapolation
end

function calculate!(curve::InterpolatedCurve)
  calculate!(curve.boot)

  return curve
end


### Fitted curve methods ###
function initialize!(curve::FittedBondDiscountCurve)
  # yield conventions
  dc = curve.dc
  yield_comp = CompoundedCompounding()
  freq = Annual()

  n = length(curve.bonds)
  cost_f = curve.fittingMethod.costFunction

  squared_sum = 0.0
  for i = 1:n
    bond = curve.bonds[i]
    leg = bond.cashflows
    clean_price = bond.faceAmount
    bond_settlement = get_settlement_date(bond)

    # get the ytm of the bond
    ytm = yield(bond, clean_price, dc, yield_comp, freq, bond_settlement)
    dur = duration(bond, ytm, dc, yield_comp, freq, ModifiedDuration(), bond_settlement)

    curve.fittingMethod.weights[i] = 1.0 / dur
    squared_sum += curve.fittingMethod.weights[i] * curve.fittingMethod.weights[i]

    cf = bond.cashflows
    for k = 1:length(cf.coupons) + 1 # for redemption
      cf_to_use = k > length(cf.coupons) ? cf.redemption : cf.coupons[i]
      if !has_occurred(cf_to_use, bond_settlement)
        cost_f.firstCashFlow[i] = k
        break
      end
    end
  end

  divide_array_by_self!(curve.fittingMethod.weights, sqrt(squared_sum))

  return curve
end

function calculate!(curve::FittedBondDiscountCurve)
  cost_f = curve.fittingMethod.costFunction
  constraint = NoConstraint()

  x = zeros(guess_size(curve.fittingMethod))

  if length(curve.fittingMethod.guessSolution) > 0
    x = curve.fittingMethod.guessSolution
  end

  simplex = Simplex(curve.simplexLambda)
  problem = Problem(cost_f, constraint, x)

  max_stationary_state_iterations = 100
  root_epsilon = curve.accuracy
  function_epsilon = curve.accuracy
  gradient_norm_epsilon = curve.accuracy

  end_criteria = EndCriteria(curve.maxEvaluations, max_stationary_state_iterations, root_epsilon, function_epsilon, gradient_norm_epsilon)

  minimize!(simplex, problem, end_criteria)
  solution = problem.currentValue

  number_of_iterations = problem.
end

function value(cf::CostFunction, x::Vector{Float64})
end
