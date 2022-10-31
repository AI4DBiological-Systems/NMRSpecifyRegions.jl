function getPsnospininfo(As::Vector{NMRHamiltonian.SHType{T}}, hz2ppmfunc) where T

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end

function evalzerophasecl1Darray(u_rad, αs::Vector{T}, Ωs::Vector{T}, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l = 1:length(αs)
        out += evalzerophaseclpartitionelement(u_rad, αs[l], Ωs[l], λ)
    end

    return out
end

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: AbstractFloat

    out = α/(λ+im*(r-Ω))

    return out
end

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

function getΩS(As::Vector{NMRHamiltonian.SHType{T}}) where T

    ΩS = Vector{Vector{Vector{T}}}(undef, length(As))

    for n = 1:length(As)

        ΩS[n] = Vector{Vector{T}}(undef, length(As[n].Ωs) + length(As[n].Ωs_singlets))
        for i = 1:length(As[n].Ωs)

            ΩS[n][i] = copy(As[n].Ωs[i])

        end

        for i = 1:length(As[n].Ωs_singlets)
            ΩS[n][length(As[n].Ωs)+i] = [ As[n].Ωs_singlets[i]; ]
        end
    end

    return ΩS
end


function getPs( ΩS::Vector{Vector{Vector{T}}}, hz2ppmfunc) where T <: Real

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

function initializeΔsyscs(As::Vector{NMRHamiltonian.SHType{T}}, x::T) where T
    N_compounds = length(As)
    Δsys_cs = Vector{Vector{T}}(undef, N_compounds)

    for n = 1:N_compounds

        N_sys = length(As[n].N_spins_sys)
        Δsys_cs[n] = Vector{T}(undef, N_sys)
        for i = 1:N_sys
            Δsys_cs[n][i] = x
        end

        for i = 1:length(As[n].αs_singlets)
            push!(Δsys_cs[n], x)
        end
    end

    return Δsys_cs
end