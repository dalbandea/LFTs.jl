
Base.@kwdef struct Phi4Parm <: LFTParm
    iL::Tuple{Int64,Int64}
    beta::Float64
    lambda::Float64
end
export Phi4Parm

@doc raw"""
    struct Phi4workspace{T}

Allocates all the necessary fields for a HMC simulation of a Phi4 model:

- `PRC`: precision; must be `Float16`, `Float32` or `Float64`.
- `phi`: ``\phi^4`` field.
- `frc`
"""
struct Phi4workspace{T, N} <: Phi4
    PRC::Type{T}
    phi::Array{T, N}
    params::Phi4Parm
    function Phi4workspace(::Type{T}, lp::Phi4Parm) where {T <: AbstractFloat}
        phi = Array{T, 2}(undef, lp.iL...)
        return new{T, 2}(T, phi, lp)
    end
end

function (::Type{Phi4})(::Type{T} = Float64; kwargs...) where {T <: AbstractFloat}
    return Phi4workspace(T, Phi4Parm(;kwargs...))
end

struct Phi4HMC{A <: AbstractArray} <: AbstractHMC
    params::HMCParams
    frc::A
    mom::A
end

function Phi4HMC(phiws::Phi4, hmcp::HMCParams)
    frc = similar(phiws.phi)
    mom = similar(phiws.phi)
    return Phi4HMC{typeof(frc)}(hmcp, frc, mom)
end

sampler(lftws::Phi4, hmcp::HMCParams) = Phi4HMC(lftws, hmcp)


function copy!(phiws_dst::Phi4, phiws_src::Phi4)
    phiws_dst.phi .= phiws_src.phi
    return nothing
end

@doc raw"""
    function randomize!(phiws, lp)

Randomizes:

- `phiws.phi` from Gaussian.
"""
function randomize!(phiws::Phi4)
    phiws.phi .= Random.randn(phiws.PRC, size(phiws.phi)...)
    return nothing
end

