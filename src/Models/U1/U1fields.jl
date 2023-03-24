
# ======================= #
# ===== General U1 ====== #
# ======================= #

function copy!(U1ws_dst::U1, U1ws_src::U1, lp::U1Parm)
    U1ws_dst.U .= U1ws_src.U
    return nothing
end

function randomize!(U1ws::U1)
    U1ws.U .= to_device(U1ws.device, exp.(im * Random.rand(U1ws.PRC, size(U1ws.U)) * 2 * pi))
    return nothing
end

# ======================= #
# ===== U1 Quenched ===== #
# ======================= #

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
    params::U1QuenchedParm
    device::Union{KernelAbstractions.Device, ROCKernels.ROCDevice}
    kprm::KernelParm
    function U1quenchedworkspace(::Type{T}, lp::U1Parm, device, kprm) where {T <: AbstractFloat}
        U = to_device(device, ones(complex(T), lp.iL..., 2))
        return new{T, typeof(U)}(T, U, lp, device, kprm)
    end
end

function (s::Type{U1Quenched})(::Type{T}; device = CUDAKernels.CUDADevice(), kwargs...) where {T <: AbstractFloat}
    lp = U1QuenchedParm(;kwargs...)
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


# ======================= #
# ====== U1 Nf = 2 ====== #
# ======================= #

struct U1Nf2workspace{T, A <: AbstractArray, S <: AbstractSolver} <: U1Nf2
    PRC::Type{T}
    U::A
    params::U1Nf2Parm
    device::Union{KernelAbstractions.Device, ROCKernels.ROCDevice}
    kprm::KernelParm
    sws::S
    function U1Nf2workspace(::Type{T}, lp::U1Nf2Parm, device, kprm, maxiter::Int64 = 10000,
            tol::Float64 = 1e-14) where {T <: AbstractFloat}
        U = to_device(device, ones(complex(T), lp.iL..., 2))
        sws = CG(maxiter, tol, U)
        return new{T, typeof(U), typeof(sws)}(T, U, lp, device, kprm, sws)
    end
end

function (s::Type{U1Nf2})(::Type{T}; device = CUDAKernels.CUDADevice(), maxiter::Int64 = 10000, tol::Float64 = 1e-14, kwargs...) where {T <: AbstractFloat}
    lp = U1Nf2Parm(;kwargs...)
    return U1Nf2workspace(T, lp, device, KernelParm(lp), maxiter, tol)
end

struct U1Nf2HMC{A1 <: AbstractArray, A2 <: AbstractArray} <: AbstractHMC
    params::HMC
    X::A1
    F::A1
    g5DX::A1
    frc1::A2 # gauge force
    frc2::A2 # gauge force
    pfrc::A2 # pf force
    mom::A2
end

function U1Nf2HMC(u1ws::U1Nf2, hmcp::HMCParams)
    X = similar(u1ws.U)
    F = similar(u1ws.U)
    g5DX = similar(u1ws.U)
    frc1 = to_device(u1ws.device, zeros(u1ws.PRC, u1ws.params.iL..., 2))
    frc2 = similar(frc1)
    pfrc = similar(frc1)
    mom = similar(frc1)
    return U1Nf2HMC{typeof(X), typeof(frc1)}(hmcp, X, F, g5DX, frc1, frc2, pfrc, mom)
end

sampler(lftws::U1Nf2, hmcp::HMCParams) = U1Nf2HMC(lftws, hmcp)


# ======================= #
# ======== U1 Nf ======== #
# ======================= #


struct U1Nfworkspace{T, A <: AbstractArray, S <: AbstractSolver} <: U1Nf
    PRC::Type{T}
    U::A
    params::U1NfParm
    device::Union{KernelAbstractions.Device, ROCKernels.ROCDevice}
    kprm::KernelParm
    sws::S
    function U1Nfworkspace(::Type{T}, lp::U1NfParm, device, kprm, maxiter::Int64 = 10000,
            tol::Float64 = 1e-14) where {T <: AbstractFloat}
        U = to_device(device, ones(complex(T), lp.iL..., 2))
        sws = CG(maxiter, tol, U)
        return new{T, typeof(U), typeof(sws)}(T, U, lp, device, kprm, sws)
    end
end
export U1Nfworkspace

function (s::Type{U1Nf})(::Type{T}; device = CUDAKernels.CUDADevice(), maxiter::Int64 = 10000, tol::Float64 = 1e-14, kwargs...) where {T <: AbstractFloat}
    lp = U1NfParm(;kwargs...)
    return U1Nfworkspace(T, lp, device, KernelParm(lp), maxiter, tol)
end

struct U1NfHMC{A1 <: AbstractArray, A2 <: AbstractArray, A3 <: AbstractArray} <: AbstractHMC
    params::HMC
    X::A1
    F::A3 # array of arrays for fermions
    g5DX::A1
    frc1::A2 # gauge force
    frc2::A2 # gauge force
    pfrc::A2 # pf force
    mom::A2
end

function U1NfHMC(u1ws::U1Nf, hmcp::HMCParams)
    X = similar(u1ws.U)
    g5DX = similar(u1ws.U)
    frc1 = to_device(u1ws.device, zeros(u1ws.PRC, u1ws.params.iL..., 2))
    frc2 = similar(frc1)
    pfrc = similar(frc1)
    mom = similar(frc1)

    N_fermions = length(u1ws.params.am0)                 # get number of fermions
    F = Array{DevArray(u1ws.device)}(undef, N_fermions)
    for i in 1:N_fermions
        F[i] = similar(u1ws.U)
    end
    return U1NfHMC{typeof(X), typeof(frc1), typeof(F)}(hmcp, X, F, g5DX, frc1, frc2, pfrc, mom)
end
export U1NfHMC

sampler(lftws::U1Nf, hmcp::HMCParams) = U1NfHMC(lftws, hmcp)
