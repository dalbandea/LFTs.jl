
abstract type QuantumRotor <: AbstractLFT end
export QuantumRotor

abstract type AbstractBoundaryCondition end
abstract type PeriodicBC <: AbstractBoundaryCondition end
abstract type AntiperiodicBC <: AbstractBoundaryCondition end
abstract type OpenBC <: AbstractBoundaryCondition end

abstract type AbstractDiscretization end
abstract type ClassicalPerfectDiscretization <: AbstractDiscretization end

struct QuantumRotorParm{B <: AbstractBoundaryCondition, D <: AbstractDiscretization} <: LFTParm
    iT::Int64
    I::Float64
    BC::Type{B}
    disc::Type{D}
end
export QuantumRotorParm

function QuantumRotorParm(; iT, I, BC::Type{B} = PeriodicBC, disc::Type{D} = ClassicalPerfectDiscretization) where {B <: AbstractBoundaryCondition, D <: AbstractDiscretization}
    return QuantumRotorParm{BC, D}(iT, I, BC, disc)
end


struct QuantumRotorWorkspace{T, N, P <: LFTParm} <: QuantumRotor
    PRC::Type{T}
    phi::Array{T, N}
    params::P
    function QuantumRotorWorkspace(::Type{T}, lp::QuantumRotorParm) where {T <: AbstractFloat}
        phi = Array{T, 1}(undef, lp.iT)
        return new{T, 1, typeof(lp)}(T, phi, lp)
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

function copy!(phiws_dst::QuantumRotor, phiws_src::QuantumRotor, hmcws::QuantumRotorHMC)
    phiws_dst.phi .= phiws_src.phi
    return nothing
end

function randomize!(phiws::QuantumRotor)
    phiws.phi .= 2pi*Random.rand(phiws.PRC, size(phiws.phi)...) .- pi
    return nothing
end


#custom mod
#changes mod domain from (0 to z) to (-z/2 to z/2)
function Mod(x, z)

    a = mod(x,z)
    if a <= z/2
        return a
    else
        return a - z
    end

end

