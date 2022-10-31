# position units are in ppm.
mutable struct MetaboliteRegionType{T}
    interval_starts::Vector{Vector{T}} # [spin_group_index][interval_index]
    interval_fins::Vector{Vector{T}}
end

# for a given compound, spin group, interval.
mutable struct RegionInfoType{T}

    # used for diagnostics.
    interval_indices::Vector{Tuple{Int,Int,Int}} # compound, spin group, interval indexing.

    # used for deciding the cost function positions for this region.
    st::T
    fin::T

    # used for updating results.
    spin_group_indices::Vector{Tuple{Int,Int}} # compound, spin group
end

mutable struct ExperimentInfoType{T}

    regions::Vector{RegionInfoType{T}}
    graph::Graphs.SimpleGraphs.SimpleGraph{Int}
    connectivity::Vector{Vector{Int}}

    sts_set::Vector{Vector{Vector{T}}}
    fins_set::Vector{Vector{Vector{T}}}
    IDs_set::Vector{Vector{Vector{Tuple{Int,Int,Int}}}}

    sts::Vector{T}
    fins::Vector{T}
    IDs::Vector{Tuple{Int,Int,Int}}
end
