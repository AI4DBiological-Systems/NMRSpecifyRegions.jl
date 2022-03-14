
function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end

"""
Returns a deepcopy of the Ωs components for each compound, and singlets.
"""
function getΩS(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}) where {T,SST}

    N = length(As)
    ΩS = Vector{Vector{Vector{T}}}(undef, N)

    for n = 1:N
        ΩS[n] = deepcopy(As[n].Ωs)

        for i = 1:length(As[n].Ωs_singlets)
            push!(ΩS[n], [As[n].Ωs_singlets[i];])
        end
    end

    return ΩS
end

function getPs( ΩS::Vector{Vector{Vector{T}}},
                hz2ppmfunc) where T <: Real

    N_compounds = length(ΩS)

    Ps = Vector{Vector{Vector{T}}}(undef, N_compounds)
    for n = 1:N_compounds

        Ps[n] = Vector{Vector{T}}(undef, length(ΩS[n]))
        for i = 1:length(ΩS[n])

            Ps[n][i] = Vector{T}(undef, length(ΩS[n][i]))
            for l = 1:length(ΩS[n][i])

                Ps[n][i][l] = hz2ppmfunc( ΩS[n][i][l]/(2*π) )
            end
        end
    end

    return Ps
end
