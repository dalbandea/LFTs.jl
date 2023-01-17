

@doc raw"""
    struct U1workspace{T}

Allocates all the necessary fields for a HMC simulation of a U(1) model:

- `PRC`: precision; must be `ComplexF64`.
- `U`: ``U`` gauge field.
- `frc`
"""
struct U1quenchedworkspace{T} <: U1Quenched
    PRC::Type{T}
    U
    frc1
    frc2
    mom
    function U1quenchedworkspace(::Type{T}, lp::U1Parm) where {T <: AbstractFloat}
        U = to_device(lp.device, ones(complex(T), lp.iL..., 2))
        frc1 = to_device(lp.device, zeros(T, lp.iL..., 2))
        frc2 = similar(frc1)
        mom = similar(frc1)
        return new{T}(T, U, frc1, frc2, mom)
    end
end

struct U1Nf2workspace{T} <: U1Nf2
    PRC::Type{T}
    U
    am0
    X
    F
    g5DX
    frc1 # gauge force
    frc2 # gauge force
    pfrc # pf force
    mom
    sws::AbstractSolver
    function U1Nf2workspace(::Type{T}, lp::U1Parm, am0, maxiter::Int64 = 10000,
            tol::Float64 = 1e-14) where {T <: AbstractFloat}
        U = to_device(lp.device, ones(complex(T), lp.iL..., 2))
        X = similar(U)
        F = similar(U)
        g5DX = similar(U)
        frc1 = to_device(lp.device, zeros(T, lp.iL..., 2))
        frc2 = similar(frc1)
        pfrc = similar(frc1)
        mom = similar(frc1)
        sws = CG(maxiter, tol, X)
        return new{T}(T, U, am0, X, F, g5DX, frc1, frc2, pfrc, mom, sws)
    end
end

function copy!(U1ws_dst::U1, U1ws_src::U1, lp::U1Parm)
    U1ws_dst.U .= U1ws_src.U
    return nothing
end

# @doc raw"""
#     function randomize!(phiws, lp)

# Randomizes:

# - `phiws.phi` from Gaussian.
# """
# function randomize!(phiws::Phi4, lp::Phi4Parm)
#     phiws.phi .= Random.randn(phiws.PRC, size(phiws.phi)...)
#     return nothing
# end

