
function getmetaboliteintervals(tag::String,
    molecule_names,
    ΩS_ppm::Vector{Vector{Vector{T}}},
    Δsys_border::Vector{Vector{T}};
    min_dist::T = 0.15) where T <: Real

    n_select = findfirst(xx->(xx==tag), molecule_names)
    C0 = ΩS_ppm[n_select]
    N_groups = length(C0)

    Δcompound_border = Δsys_border[n_select]
    @assert N_groups == length(Δcompound_border)

    st_c = Vector{Vector{T}}(undef, N_groups)
    fin_c = Vector{Vector{T}}(undef, N_groups)

    for i = 1:N_groups
        st_c[i], fin_c[i] = segmentmoleculefreqs2(C0[i],
            Δcompound_border[i], min_dist = min_dist)
    end

    return MetaboliteRegionType(st_c, fin_c)
end

# A region needs to contain all interval C0[l] +/- r such that
#   all intervals are within min_dist from each other.
# r and min_dst are in ppm.
function segmentmoleculefreqs2(C0::Vector{T}, border::T;
    min_dist::T = 0.15) where T <: Real

    sort_inds = sortperm(C0)
    C = C0[sort_inds]

    st_c = Vector{T}(undef, 0)
    fin_c = Vector{T}(undef, 0)

    push!(st_c, C[1])

    for l = 2:length(C)
        left = C[l-1] + border
        right = C[l] -border

        if (right-left) > min_dist
            push!(fin_c, C[l-1])
            push!(st_c, C[l])
        end
    end
    push!(fin_c, C[end])

    return st_c, fin_c
end

function constructgraph(sts::Vector{T}, fins::Vector{T}) where T

    # construct graph.
    N_intervals = length(sts)
    g = LightGraphs.SimpleGraph(N_intervals)

    for i = 1:N_intervals
        flags = isthereoverlap(sts[i], fins[i], sts, fins)

        for j = 1:length(flags)
            if flags[j]
                #add edge i to j, j to i.
                LightGraphs.add_edge!(g, i, j)
                LightGraphs.add_edge!(g, j, i)
            end
        end
    end

    return g
end

function isthereoverlap(st::T, fin::T, sts, fins) where T

    flags = collect( isoverlapped(st, fin, sts[i], fins[i]) for i = 1:length(fins))

    return flags
end


# assumes st_a < fin_a, st_b < fin_b.
function isoverlapped(st_a::T, fin_a::T, st_b::T, fin_b::T) where T

    if st_a < st_b && st_b < fin_a < fin_b
        return true

    elseif st_b < st_a && st_a < fin_b < fin_a
        return true

    elseif st_b < st_a < fin_b
        return true

    elseif st_a < st_b < fin_a
        return true

    end

    return false
end
