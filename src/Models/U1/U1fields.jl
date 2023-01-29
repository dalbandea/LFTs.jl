

@doc raw"""
    struct U1workspace{T}

Allocates all the necessary fields for a HMC simulation of a U(1) model:

- `PRC`: precision; must be `ComplexF64`.
- `U`: ``U`` gauge field.
- `frc`
"""
struct U1quenchedworkspace{T, A <: AbstractArray} <: U1Quenched
    PRC::Type{T}
    U::A
    params::U1Parm
    device::Union{KernelAbstractions.Device, ROCKernels.ROCDevice}
    kprm::KernelParm
    function U1quenchedworkspace(::Type{T}, lp::U1Parm, device, kprm) where {T <: AbstractFloat}
        U = to_device(device, ones(complex(T), lp.iL..., 2))
        return new{T, typeof(U)}(T, U, lp, device, kprm)
    end
end

function (s::Type{U1Quenched})(::Type{T}; device = CUDAKernels.CUDADevice(), kwargs...) where {T <: AbstractFloat}
    lp = U1Parm(;kwargs...)
    return U1quenchedworkspace(T, lp, device, KernelParm(lp))
end

struct U1quenchedHMC{A <: AbstractArray} <: AbstractHMC
    params::HMC
    frc1::A
    frc2::A
    mom::A
end


function U1quenchedHMC(u1ws::U1Quenched, hmcp::HMCParams)
    frc1 = to_device(u1ws.device, zeros(u1ws.PRC, u1ws.params.iL..., 2))
    frc2 = similar(frc1)
    mom = similar(frc1)
    return U1quenchedHMC(hmcp, frc1, frc2, mom)
end

sampler(lftws::U1Quenched, hmcp::HMCParams) = U1quenchedHMC(lftws, hmcp)

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

function randomize!(U1ws::U1, lp::U1Parm)
    U1ws.U .= exp.(im*(2pi*Random.rand(U1ws.PRC, size(U1ws.U)...) .- pi))
    return nothing
end

