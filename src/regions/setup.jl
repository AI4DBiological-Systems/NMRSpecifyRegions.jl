function setupexperimentresults(molecule_names,
    ΩS_ppm,
    Δsys_border::Vector{Vector{T}};
    min_dist = 0.15) where T

    metabolite_regions = collect( getmetaboliteintervals(molecule_names[n],
        molecule_names, ΩS_ppm, Δsys_border;
        min_dist = min_dist) for n = 1:length(molecule_names) )

    # get connected graph. A region is defined as a connected subgraph.
    D, sts_set, fins_set, IDs_set,
        G = setupregioninfo(metabolite_regions, Δsys_border)
    # D is the list of graph nodes that form connected components.

    # information from each graph nodes.
    sts = combinevectors(combinevectors(sts_set))
    fins = combinevectors(combinevectors(fins_set))
    IDs = combinevectors(combinevectors(IDs_set))

    regions = Vector{RegionInfoType{T}}(undef, length(D))
    for i = 1:length(D)

        inds = D[i]
        ids = IDs[inds]

        s = collect( ids[k][1:2] for k = 1:length(ids) )
        spin_group_indices = sort(unique(s))

        # find start and fin of region.
        tmp = collect( fins[inds[k]] for k = 1:length(inds) )
        region_fin = maximum(tmp)

        tmp = collect( sts[inds[k]] for k = 1:length(inds) )
        region_st = minimum(tmp)

        regions[i] = RegionInfoType(ids, region_st, region_fin, spin_group_indices)
    end

    return ExperimentInfoType(regions, G, D, sts_set, fins_set, IDs_set,
        sts, fins, IDs)
end

function setupregioninfo(metabolite_regions,
    Δsys_border::Vector{Vector{T}}) where T

    # get connected graph.
    sts_set, fins_set, IDs_set, G = processintervals(metabolite_regions, Δsys_border)
    D0 = LightGraphs.connected_components(G)

    sts = combinevectors(combinevectors(sts_set))
    fins = combinevectors(combinevectors(fins_set))
    IDs = combinevectors(combinevectors(IDs_set))

    # sort D by size.
    region_sizes = collect( length(D0[i]) for i = 1:length(D0))
    inds = sortperm(region_sizes)
    region_sizes = region_sizes[inds]
    D = D0[inds]

    # sanity-check.
    @assert length(unique(combinevectors(D))) == LightGraphs.nv(G)

    # schedule jobs.
    IDs = combinevectors(combinevectors(IDs_set))

    return D, sts_set, fins_set, IDs_set, G
end


# border in ppm, adds that amount to each metabolite interval.
function processintervals(metabolite_regions::Vector{MetaboliteRegionType{T}},
    Δsys_border::Vector{Vector{T}}) where T


    # initialize regions to every interval.
    N_compounds = length(metabolite_regions)
    @assert length(Δsys_border) == N_compounds

    sts_set = Vector{Vector{Vector{T}}}(undef, N_compounds)
    fins_set = Vector{Vector{Vector{T}}}(undef, N_compounds)
    IDs_set = Vector{Vector{Vector{Tuple{Int,Int,Int}}}}(undef, N_compounds)

    for n = 1:N_compounds

        N_groups = length(Δsys_border[n])

        sts_set[n] = Vector{Vector{T}}(undef, N_groups)
        fins_set[n] = Vector{Vector{T}}(undef, N_groups)
        IDs_set[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_groups)

        for i = 1:N_groups

            sts_set[n][i] = metabolite_regions[n].interval_starts[i] .- Δsys_border[n][i]
            fins_set[n][i] = metabolite_regions[n].interval_fins[i] .+ Δsys_border[n][i]

            N_intervals = length( metabolite_regions[n].interval_starts[i] )
            #IDs_set[n][i] = (n,i,length(metabolite_regions[n].interval_fins[i]))
            IDs_set[n][i] = collect( (n,i,k) for k = 1:N_intervals )
        end
    end

    sts = combinevectors(combinevectors(sts_set))
    fins = combinevectors(combinevectors(fins_set))

    # sanity check.
    M2 = getnumelements(sts_set)
    M3 = getnumelements(IDs_set)
    @assert M2 == M3

    # construct graph.
    G = constructgraph(sts, fins)

    return sts_set, fins_set, IDs_set, G
end

# get total number of elements.
function getnumelements(X::Vector{Vector{Vector{T}}})::Int where T
    return sum( sum( length(X[n][i]) for i = 1:length(X[n]) ) for n = 1:length(X) )
end


#########

function filterfreqpositions(U_cost,
    st_us::Vector{T}, fin_us::Vector{T}) where T

    band_inds = Vector{Int}(undef, 0)
    for m = 1:length(st_us)

        st_u = st_us[m]
        fin_u = fin_us[m]
        inds = findall(xx->(st_u<xx<fin_u), U_cost)
        push!(band_inds, inds...)
    end

    return unique(band_inds)
end
