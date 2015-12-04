# Solvers

abstract Solver1D

type BrentSolver <: Solver1D end

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

# brent solver function
function _solve(solver::BrentSolver, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64, root::Float64, eval_num::Int64)
  froot = f(root)
  eval_num += 1
  max_evals = 100 # const?
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
