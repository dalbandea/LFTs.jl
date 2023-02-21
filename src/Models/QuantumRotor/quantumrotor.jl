
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

include("quantumrotorfields.jl")
include("quantumrotoraction.jl")
include("quantumrotorhmc.jl")
include("quantumrotormeasurements.jl")
