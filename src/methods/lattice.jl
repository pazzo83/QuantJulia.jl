type TreeLattice1D{I <: Integer} <: TreeLattice
  tg::TimeGrid
  statePrices::Vector{Vector{Float64}}
  n::I
  statePricesLimit::I
end

function TreeLattice1D{I}(tg::TimeGrid, n::I)
  statePrices = Vector{Vector{Float64}}(1)
  statePrices[1] = ones(1)

  statePricesLimit = 0

  return TreeLattice1D(tg, statePrices, n, statePricesLimit)
end

type Branching{I <: Integer}
  k::Vector{I}
  probs::Vector{Vector{Float64}}
  kMin::I
  jMin::I
  kMax::I
  jMax::I
end

Branching() = Branching(zeros(Int, 0), Vector{Vector{Float64}}(3), typemax(Int), typemax(Int), typemin(Int), typemin(Int))

get_size(b::Branching) = b.jMax - b.jMin + 1

function add!(branch::Branching, k::Float64, p1::Float64, p2::Float64, p3::Float64)
  push!(branch.k, k)
  push!(branch.probs[1], p1)
  push!(branch.probs[2], p2)
  push!(branch.probs[3], p3)

  # maintain invariants
  branch.kMin = min(branch.kMin, k)
  branch.jMin = branch.kMin - 1
  branch.kMax = max(branch.kMax, k)
  branch.jMax = branch.kMax - 1

  return branch
end

type TrinomialTree{S <: StochasticProcess}
  process::S
  timeGrid::TimeGrid
  dx::Vector{Float64}
  branchings::Vector{Branching}
  isPositive::Bool

  function TrinomialTree{S}(process::S, timeGrid::TimeGrid, isPositive::Bool = false)
    x0 = process.x0
    dx = zeros(length(timeGrid.times))
    nTimeSteps = length(timeGrid.times) - 1
    jMin = 1
    jMax = 1
    branchings = Vector{Branching}(nTimeSteps)

    for i = 1:nTimeSteps
      t = timeGrid.times[i]
      dt = timeGrid.dt[i]

      # Variance must be independent of x
      v2 = variance(process, t, 0.0, dt)
      v = sqrt(v2)
      dx[i+1] = v * sqrt(3.0)

      branching = Branching()

      for i =jMin:jMax
        x = x0 + j * dx[i]
        m = expectation(process, t, x, dt)
        temp = round(Int, floor((m - x0) / dx[i+1] + 0.5))

        if isPositive
          while (x0 + (temp - 1) * dx[i + 1] <= 0)
            temp += 1
          end
        end

        e = m - (x0 + temp * dx[i + 1])
        e2 = e * e
        e3 = e * sqrt(3.0)

        p1 = (1.0 + e2 / v2 - e3 / v) / 6.0
        p2 = (2.0 - e2 / v2) / 3.0
        p3 = (1.0 + e2 / v2 + e3 / v) / 6.0

        add!(branching, temp, p1, p2, p3)
      end

      branchings[i] = copy(branching) # check if we need copy

      jMin = branching.jMin
      jMax = branchin.jMax
    end
  end
end

get_size{I <: Integer}(t::TrinomialTree, n::I) = i == 0 ? 1 : get_size(t.branchings[i])
