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

type FittedBondDiscountCurve{B <: Bond} <: Curve
  settlementDays::Integer
  referenceDate::Date
  calendar::BusinessCalendar
  bonds::Vector{B}
  dc::DayCount
  fittingMethod::FittingMethod
  accuracy::Float64
  maxEvaluations::Integer
  simplexLambda::Float64

  FittedBondDiscountCurve(settlementDays::Integer,
                          referenceDate::Date,
                          calendar::BusinessCalendar,
                          bonds::Vector{B},
                          dc::DayCount,
                          fittingMethod::FittingMethod,
                          accuracy::Float64,
                          maxEvaluations::Integer,
                          simplexLambda::Float64) =

                          (x = new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda);
                          x.fittingMethod.commons.costFunction.curve = x)

  # FittedBondDiscountCurve(settlementDays::Integer, referenceDate::Date, calendar::BusinessCalendar, bonds::Vector{B}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64=1e-10,
  #                                     maxEvaluations::Integer=10000, simplexLambda::Float64=1.0) =
  #                         new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
    #n = length(bonds)
  #   println("hi")
  #   this = new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
  #   println("also hi")
  #   this.fittingMethod.costFunction.curve = this
  #
  #   return this
  # end
end

FittedBondDiscountCurve{B <: Bond}(settlementDays::Integer, referenceDate::Date, calendar::BusinessCalendar, bonds::Vector{B}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64=1e-10,
                                     maxEvaluations::Integer=10000, simplexLambda::Float64=1.0) = FittedBondDiscountCurve{B}(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)

type FittingCost <: CostFunction
  # value::Vector{Float64}
  # values::Vector{Float64}
  firstCashFlow::Vector{Integer}
  curve::Curve
end

function FittingCost(size::Integer, curve::Curve)
  # value = Vector{Float64}()
  # values = Vector{Float64}()
  firstCashFlow = zeros(Int, size)

  return FittingCost(firstCashFlow, curve)
end

# Interpolated Curve methods #
max_date(curve::InterpolatedCurve) = curve.dates[end]

function discount(curve::Curve, t::Float64)
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
discount_impl(curve::FittedBondDiscountCurve, t::Float64) = discount_function(curve.fittingMethod, curve.fittingMethod.commons.solution, t)

function initialize!(curve::FittedBondDiscountCurve)
  # yield conventions
  dc = curve.dc
  yield_comp = CompoundedCompounding()
  freq = Annual()

  n = length(curve.bonds)
  cost_f = curve.fittingMethod.commons.costFunction

  squared_sum = 0.0
  for i = 1:n
    bond = curve.bonds[i]
    leg = bond.cashflows
    clean_price = bond.faceAmount
    bond_settlement = get_settlement_date(bond)

    # get the ytm of the bond
    ytm = yield(bond, clean_price, dc, yield_comp, freq, bond_settlement)
    dur = duration(bond, ytm, dc, yield_comp, freq, ModifiedDuration(), bond_settlement)

    curve.fittingMethod.commons.weights[i] = 1.0 / big(dur)
    squared_sum += curve.fittingMethod.commons.weights[i] * curve.fittingMethod.commons.weights[i]

    cf = bond.cashflows
    for k = 1:length(cf.coupons) + 1 # for redemption
      cf_to_use = k > length(cf.coupons) ? cf.redemption : cf.coupons[i]
      if !has_occurred(cf_to_use, bond_settlement)
        cost_f.firstCashFlow[i] = k
        break
      end
    end
  end

  divide_array_by_self!(curve.fittingMethod.commons.weights, sqrt(squared_sum))

  return curve
end

function calculate!(curve::FittedBondDiscountCurve)
  cost_f = curve.fittingMethod.commons.costFunction
  constraint = NoConstraint()

  x = zeros(curve.fittingMethod.size)

  if length(curve.fittingMethod.commons.guessSolution) > 0
    x = curve.fittingMethod.commons.guessSolution
  end

  simplex = Simplex(curve.simplexLambda)
  problem = Problem(cost_f, constraint, x)

  max_stationary_state_iterations = 100
  root_epsilon = curve.accuracy
  function_epsilon = curve.accuracy
  gradient_norm_epsilon = curve.accuracy

  end_criteria = EndCriteria(curve.maxEvaluations, max_stationary_state_iterations, root_epsilon, function_epsilon, gradient_norm_epsilon)

  minimize!(simplex, problem, end_criteria)
  curve.fittingMethod.commons.solution = problem.currentValue

  number_of_iterations = problem.functionEvaluation
  cost_value = problem.functionValue

  curve.fittingMethod.commons.guessSolution = curve.fittingMethod.commons.solution
  curve.fittingMethod.commons.numberOfIterations = number_of_iterations
  curve.fittingMethod.commons.minimumCostValue = cost_value

  return curve
end

function value{T}(cf::CostFunction, x::Vector{T})
  ref_date = cf.curve.referenceDate
  dc = cf.curve.dc
  squared_error = 0.0
  n = length(cf.curve.bonds)

  for (i, bond) in enumerate(cf.curve.bonds)
    bond_settlement = get_settlement_date(bond)
    model_price = -accrued_amount(bond, bond_settlement)
    leg = bond.cashflows
    for k = cf.firstCashFlow[i]:length(leg.coupons)
      # @inbounds df = discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
      @inbounds model_price += amount(leg.coupons[k]) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
    end

    # redemption
    @inbounds model_price += amount(leg.redemption) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.redemption)))

    # adjust NPV for forward settlement
    if bond_settlement != ref_date
      model_price /= discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, bond_settlement))
    end

    market_price = bond.faceAmount
    price_error = model_price - market_price
    @inbounds weighted_error = cf.curve.fittingMethod.commons.weights[i] * price_error
    squared_error += weighted_error * weighted_error
  end

  return squared_error
end
