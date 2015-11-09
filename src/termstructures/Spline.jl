type Spline
  x_vert::Vector{Int64}
  y_vert::Vector{Float64}
  b::Vector{Float64}
  c::Vector{Float64}
  d::Vector{Float64}
end

function gen_splines(x_vert::Vector{Int64}, y_vert::Vector{Float64})
  n = length(x_vert)
  h = zeros(n - 1)
  b = zeros(n - 1)
  a = zeros(n - 2)
  g = zeros(n - 2)
  c = zeros(n)
  d = zeros(n - 1)

  for i = 1:n - 1
    h[i] = x_vert[i + 1] - x_vert[i]
    b[i] = (y_vert[i + 1] - y_vert[i]) / h[i]
  end

  for i = 1:n - 2
    a[i] = 2.0 * (h[i] + h[i + 1])
    # u[i] = 6 * (b[i + 1] - b[i])
    g[i] = b[i + 1] - b[i]
  end

  Alu = lufact(Tridiagonal(h[2:end-1], a[1:end], h[2:end-1]))
  c[2:end - 1] = Alu \ g

  for i=1:n - 1
    d[i] = (c[i+1] - c[i])/h[i]
    b[i] -= (2.0 * c[i] + c[i + 1]) * h[i]
    c[i] *= 3.0
  end

  return Spline(x_vert, y_vert, b, c, d)
end

function spl_interp(spl::Spline, my_x::Integer)
  # get coefficients
  # b, c, d = gen_splines(x_vert, y_vert)

  x_idx = searchsortedlast(spl.x_vert, my_x)

  if x_idx == 0
    x_idx = 1
  elseif x_idx == length(spl.x_vert)
    x_idx = length(spl.x_vert) - 1
  end


  diff = my_x - spl.x_vert[x_idx]

  return spl.y_vert[x_idx] + diff*(spl.b[x_idx] + diff * (spl.c[x_idx] + (diff * spl.d[x_idx])))
end
