# QuantJulia

module QuantJulia

# Time module
include("time/Time.jl")

# Math module
include("math/Math.jl")

# Quotes
include("quotes/Quotes.jl")

# Interest Rates
include("InterestRates.jl")

# Term Structures
include("termstructures/TermStructures.jl")

# cash flows
include("cash_flows/cash_flows.jl")

# Instruments
include("instruments/instruments.jl")

# Pricing Engines
include("pricing_engines/pricing_engines.jl")

# # Helpers NOW IN TERM STRUCTURE
# include("helpers/bond_helpers.jl")

type Settings
  evaluation_date::Date
end

settings = Settings(Date())

function set_eval_date!(sett::Settings, d::Date)
  sett.evaluation_date = d
end

export Settings, settings, set_eval_date!

using QuantJulia.Time, QuantJulia.Quotes, QuantJulia.TermStructures

end
