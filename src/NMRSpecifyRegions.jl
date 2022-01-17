module NMRSpecifyRegions

# Write your package code here.
import NMRSpectraSimulator
import LightGraphs

include("../src/types.jl")
include("../src/utils.jl")

include("../src/IO/string_helpers.jl")
include("../src/IO/config_helpers.jl")


include("../src/regions/setup.jl")
include("../src/regions/overlap.jl")

end
