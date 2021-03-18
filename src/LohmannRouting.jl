__precompile__()

module LohmannRouting

include("utils.jl")
include("watersheds.jl")
include("routing.jl")

export upstream_catchments,
       flow_sequence
       calc_xmask,
       findnearest,
       upsample,
       match_dates,
       make_irf,
       make_river_irf,
       make_grid_uh,
       rout

end # module
