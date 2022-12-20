
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
    frc::Array{T, N}
    mom::Array{T, N}
    function Phi4workspace(::Type{T}, lp::Phi4Parm) where {T <: AbstractFloat}
        phi = Array{T, 2}(undef, lp.iL...)
        frc = similar(phi)
        mom = similar(phi)
        return new{T, 2}(T, phi, frc, mom)
    end
end

function copy!(phiws_dst::Phi4, phiws_src::Phi4, lp::Phi4Parm)
    phiws_dst.phi .= phiws_src.phi
    return nothing
end

@doc raw"""
    function randomize!(phiws, lp)

Randomizes:

- `phiws.phi` from Gaussian.
"""
function randomize!(phiws::Phi4, lp::Phi4Parm)
    phiws.phi .= Random.randn(phiws.PRC, size(phiws.phi)...)
    return nothing
end

