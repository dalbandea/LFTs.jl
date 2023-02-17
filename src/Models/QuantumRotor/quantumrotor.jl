
abstract type QuantumRotor <: AbstractLFT end
export QuantumRotor

abstract type AbstractBoundaryCondition end
abstract type PeriodicBC <: AbstractBoundaryCondition end
abstract type AntiperiodicBC <: AbstractBoundaryCondition end
abstract type OpenBC <: AbstractBoundaryCondition end

abstract type AbstractDiscretization end

Base.@kwdef struct QuantumRotorParm{B <: AbstractBoundaryCondition} <: LFTParm
    iT::Int64
    I::Float64
    BC::B
end
export QuantumRotorParm

function QuantumRotorParm(; iT, I, BC::Type{B} = PeriodicBC) where {B <: AbstractBoundaryCondition}
    return QuantumRotorParm(iT = iT, I = I, BC = BC)
end


struct QuantumRotorWorkspace{T, N} <: QuantumRotor
    PRC::Type{T}
    phi::Array{T, N}
    params::Phi4Parm
    function QuantumRotorWorkspace(::Type{T}, lp::QuantumRotorParm) where {T <: AbstractFloat}
        phi = Array{T, 1}(undef, lp.iT)
        return new{T, 1}(T, phi, lp)
    end
end

function (::Type{QuantumRotor})(::Type{T} = Float64; kwargs...) where {T <: AbstractFloat}
    return QuantumRotorWorkspace(T, QuantumRotorParm(;kwargs...))
end



struct QuantumRotorHMC{A <: AbstractArray} <: AbstractHMC
    params::HMCParams
    frc::A
    mom::A
end

function QuantumRotorHMC(phiws::QuantumRotor, hmcp::HMCParams)
    frc = similar(phiws.phi)
    mom = similar(phiws.phi)
    return QuantumRotorHMC{typeof(frc)}(hmcp, frc, mom)
end

sampler(lftws::QuantumRotor, hmcp::HMCParams) = QuantumRotorHMC(lftws, hmcp)

function copy!(phiws_dst::QuantumRotor, phiws_src::QuantumRotor)
    phiws_dst.phi .= phiws_src.phi
    return nothing
end

function randomize!(phiws::QuantumRotor)
    phiws.phi .= Random.randn(phiws.PRC, size(phiws.phi)...)
    return nothing
end
