# Solvers

abstract Solver1D

type SolverInfo
  maxEvals::Integer
  lowerBoundEnforced::Bool
  upperBoundEnforced::Bool
  lowerBound::Float64
  upperBound::Float64
end


type BrentSolver <: Solver1D
  solverInfo::SolverInfo
end
BrentSolver(maxEvals::Integer = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0) =
  BrentSolver(SolverInfo(maxEvals, lowerBoundEnforced, upperBoundEnforced, lowerBound, upperBound))

type NewtonSolver <: Solver1D
  solverInfo::SolverInfo
end
NewtonSolver(maxEvals::Integer = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0) =
  NewtonSolver(SolverInfo(maxEvals, lowerBoundEnforced, upperBoundEnforced, lowerBound, upperBound))

type FiniteDifferenceNewtonSafe <: Solver1D
  solverInfo::solverInfo
end
FiniteDifferenceNewtonSafe(maxEvals::Integer = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0) =
  FiniteDifferenceNewtonSafe(solverInfo(maxEvals, lowerBoundEnforced, upperBoundEnforced, lowerBound, upperBound))

# misc functions for solving
function is_close(x::Float64, y::Float64, n::Int64 = 42)
  if x == y
    return true
  end

  diff = abs(x - y)
  tol = n * eps(Float64) # machine epsilon

  if (x * y == 0.0) # x or y is 0
    return diff < (tol * tol)
  end

  return diff <= tol * abs(x) && diff <= tol * abs(y)
end

# solver functions
function solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, step::Float64)
  ## This method returns the 0 of a function determined by a given accuracy.  This method using bracketing
  ## routine to which an intial guess must be supplied as well as a step
  growth_factor = 1.6
  flipflop = -1

  root = guess
  fxMax = f(root)

  # monotonically crescent bias, as in optionValue(volatility)
  if is_close(fxMax, 0.0)
    return root
  else if fxMax > 0.0
    xMin = enforced_bounds(solver, root - step)
    fxMin = f(xMin)
    xMax = root
  else
    xMin = root
    fxMin = fxMax
    xMax = enforced_bounds(solver, root + step)
    fxMax = f(xMax)
  end

  eval_num = 2

  while eval_num < solver.solverInfo.maxEvals
    if fxMin * fxMax <= 0.0
      if is_close(fxMin, 0.0)
        return xMin
      end
      if is_close(fxMax, 0.0)
        return xMax
      end
      root = (xMax + xMin) / 2.0
      return _solve(solver, f, accuracy, xMin, xMax, root, eval_num)
    end

    if abs(fxMin) < abs(fxMax)
      xMin = enforced_bounds(solver, xMin + growth_factor * (xMin - xMax))
      fxMin = f(xMin)
    else if abs(fxMin) > abs(fxMax)
      xMax = enforced_bounds(solver, xMax + growth_factor * (xMax - xMin))
      fxMax = f(xMax)
    else if flipflop == -1
      xMin = enforced_bounds(solver, xMin + growth_factor * (xMin - xMax))
      fxMin = f(xMin)
      eval_num += 1
      flipflop = 1
    else if flipflop == 1
      xMax = enforced_bounds(solver, xMax + growth_factor * (xMax - xMin))
      fxMax = f(xMax)
      flipflop = -1
    end

    eval_num += 1
  end

  error("Cannot converge!")
end

function solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, xMin::Float64, xMax::Float64)
  fxMin = f(xMin)
  if is_close(fxMin, 0.0)
    return xMin
  end

  fxMax = f(xMax)
  if is_close(fxMax, 0.0)
    return xMax
  end

  eval_num = 2

  return _solve(solver, f, accuracy, xMin, xMax, fxMin, fxMax, guess, eval_num)
end

### INDIVIDUAL SOLVER IMPLEMENTATIONS ###
# brent solver function
function _solve(solver::BrentSolver, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64, root::Float64, eval_num::Integer)
  froot = f(root)
  eval_num += 1
  max_evals = solver.solverInfo.maxEvals
  if froot * fxMin < 0
    xMax = xMin
    fxMax = fxMin
  else
    xMin = xMax
    fxMin = fxMax
  end

  d = root - xMax
  e = d

  while eval_num < max_evals
    if (froot > 0.0 && fxMax > 0.0) || (froot < 0.0 && fxMax < 0.0)
      # rename xMin, root, xMax, and adjust the bounds
      xMax = xMin
      fxMax = fxMin
      e = d = root - xMin
    end

    if abs(fxMax) < abs(froot)
      xMin = root
      root = xMax
      xMax = xMin
      fxMin = froot
      froot = fxMax
      fxMax = fxMin
    end

    # convergence check
    xAcc1 = 2.0 * eps(Float64) * abs(root) + 0.5 * accuracy
    xMid = (xMax - root) / 2.0
    if abs(xMid) <= xAcc1 || is_close(froot, 0.0)
      root = f(root)
      eval_num += 1
      return root
    end

    if abs(e) >= xAcc1 && abs(fxMin) > abs(froot)
      # attempting inverse quadratic interpolation
      s = froot / fxMin
      if is_close(xMin, xMax)
        p = 2.0 * xMid * s
        q = 1.0 - s
      else
        q = fxMin / fxMax
        r = froot / fxMax
        p = s * (2.0 * xMid * q * (q - r) - (root - xMin) * (r - 1.0))
        q = (q - 1.0) * (r - 1.0) * (s - 1.0)
      end

      if p > 0.0
        q = -q # check if in bounds
      end

      p = abs(p)
      min1 = 3.0 * xMid * q - abs(xAcc1 * q)
      min2 = abs(e * q)

      if (2.0 * p < (min1 < min2 ? min1 : min2))
        e = d
        d = p / q
      else
        # interpolation failed, used bisection
        d = xMid
        e = d
      end
    else
      # bounds are increasing too slowly
      d = xMid
      e = d
    end
    xMin = root
    fxMin = froot
    if abs(d) > xAcc1
      root += d
    else
      root += sign(xMid) * xAcc1
    end

    froot = f(root)
    eval_num += 1
  end

  error("Maximum number of function evals exceeded!")
end

# Finite Differences Solver
function _solve(solver::FiniteDifferenceNewtonSafe, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64,
                root::Float64, eval_num::Integer)

  max_evals = solver.solverInfo.maxEvals
  # orienting search such that f(xl) < 0
  if fxMin < 0.0
    xl = xMin
    xh = xMax
  else
    xh = xMin
    xl = xMax
  end

  froot = f(root)
  eval_num += 1

  # first order finite difference derivative
  dfroot = xMax - root < root - xMin ? (fxMax - froot) / (xMax - root) : (fxMin - froot) / (xMin - root)

  dx = xMax - xMin
  while eval_num <= max_evals
    frootold = froot
    rootold = root
    dxold = dx
    # bisect if out of range or not decreasing fast enough
    if (((root - xh) * dfroot - froot) * ((root - xl) * dfroot - froot)) > 0.0 || abs(2.0 * froot) > abs(dxold * dfroot)
      dx = (xh - xl) / 2.0
      root = xl + dx
      # if root estimate computed is close to previous one, calculate dfroot at root ane xh rather than
      # root and rootold
      if is_close(root, rootold, 2500)
        rootold = xh
        frootold = f(xh)
      end
    else
      # Newton
      dx = froot / dfroot
      root -= dx
    end

    # Convergence check
    if abs(dx) < accuracy
      return root
    end

    froot = f(root)
    eval_num += 1
    dfroot = (frootold - froot) / (rootold - root)

    if (froot < 0.0)
      xl = root
    else
      xh = root
    end
  end

  error("Maximum number of function evals exceeded!")
end

# Newton Solver (safe)
function _solve(solver::NewtonSolver, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64, root::Float64, eval_num::Integer)
  # Orient the search so that f(xl) < 0
  if fxMin < 0.0
    xl = xMin
    xh = xMax
  else
    xh = xMin
    xl = xMax
  end

  # the stepsize before last
  dxold = xMax - xMin

  # and the last step
  dx = dxold

  froot = f(root)
  dfroot = f(root, Derivative()) # derivative, hurray multiple dispatch!
  eval_num += 1

  while eval_num < solver.solverInfo.maxEvals
    # bisect if out of range or not decreasing fast enough
    if ((root - xh) * dfroot - froot) * ((root - xl) * dfroot - froot) > 0.0 || abs(2.0 * froot) > abs(dxold * dfroot)
      dxold = dx
      dx = (xh - xl) / 2.0
      root = xl + dx
    else
      dxold = dx
      dx = froot / dfroot
      root -= dx
    end

    # Convergence criterion
    if abs(dx) < accuracy
      root = f(root) # check this
      eval_num += 1
      return root
    end

    froot = f(root)
    dfroot = f(root, Derivative())
    eval_num += 1
    if froot < 0.0
      xl = root
    else
      xh = root
    end
  end

  error("Maximum number of function evals exceeded!")
end

function enforced_bounds(solver::Solver1d, x::Float64)
  if solver.solverInfo.lowerBoundEnforced && x < solver.solverInfo.lowerBound
    return solver.solverInfo.lowerBound
  end

  if solver.solverInfo.upperBoundEnforced && x > solver.solverInfo.upperBound
    return solver.solverInfo.upperBound
  end

  return x
end
