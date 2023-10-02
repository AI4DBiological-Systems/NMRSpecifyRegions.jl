function extractinfofromconfig(
    ::Type{T},
    config_path::String,
    molecule_names::Vector{String},
    ) where T <: Real

    N_compounds = length(molecule_names)

    # read.
    file_strings = readlines(config_path)

    cs_delta_group = Vector{Vector{Vector{T}}}(undef, N_compounds)
    #λ_group_labels = Vector{Vector{Vector{Int}}}(undef, N_compounds)

    for n = 1:N_compounds

        cs_delta_group[n] = extractmoleculeinfofromconfig(
            T,
            file_strings,
            molecule_names[n],
            config_path,
        )
    end

    return cs_delta_group#, λ_group_labels
end

function condenseΔcsconfig(cs_delta_group::Vector{Vector{Vector{T}}}) where T
    #
    out = Vector{Vector{T}}(undef, length(cs_delta_group))

    for n = 1:length(cs_delta_group)
        out[n] = Vector{T}(undef, length(cs_delta_group[n]))

        for i = 1:length(cs_delta_group[n])
            out[n][i] = cs_delta_group[n][i][1]
        end
    end

    return out
end


function extractmoleculeinfofromconfig(
    ::Type{T},
    file_strings,
    name::String,
    file_label::String,
    ) where T <: Real

    # find the block of text for the input name.
    text_buffer = findblock(file_strings, name, "end compound", file_label)
    if isempty(text_buffer)
        println("Error processing config file.")
        return Vector{Vector{T}}(undef, 0), Vector{Vector{Int}}(undef, 0)
    end

    # get the number of groups.
    N_groups = length(filter(xx->occursin("Group",xx), text_buffer))

    # get default set up.
    cs_delta_group = Vector{Vector{T}}(undef, N_groups)
    #λ_group_labels = Vector{Vector{Int}}(undef, N_groups)

    for i = 1:N_groups

        # update position.
        cs_buffer = findblock(text_buffer, "Group $(i)", "end group", "Group $(i)")
        if isempty(cs_buffer)
            return Vector{Vector{T}}(undef, 0), Vector{Vector{Int}}(undef, 0)
        end
        N_cs = length(cs_buffer)

        cs_delta_group[i] = Vector{T}(undef, N_cs)
        #λ_group_labels[i] = Vector{Int}(undef, N_cs)

        for k in eachindex(cs_buffer)
            tokens = split(cs_buffer[k], (':',','))
            cs_delta_group[i][k] = tryparse(T, tokens[2])
            #λ_group_labels[i][k] = tryparse(Int, tokens[3])
        end
    end

    return cs_delta_group#, λ_group_labels
end
