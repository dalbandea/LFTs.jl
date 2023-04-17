abstract type AbstractSampler end
abstract type SamplerParameters end

abstract type AbstractHMC <: AbstractSampler end
abstract type HMCParams <: SamplerParameters end

abstract type AbstractHMCLFT end

function sampler(lftws::AbstractLFT, samplerparms::SamplerParameters) end
function sample!(lftws::AbstractLFT, samplerws::AbstractSampler) end 
function sample!(lftws::AbstractLFT, samplerparms::SamplerParameters)
    samplerws = sampler(lftws, samplerparms)
    return sample!(lftws, samplerws)
end
export sample!

include("HMC/integrators/integrators.jl")
export Leapfrog, OMF4
include("HMC/hmc.jl")
export hmc!

Base.@kwdef mutable struct HMC <: HMCParams
    integrator::AbstractIntegrator = Leapfrog()
    ntherm::Int64 = 10
    ntraj::Int64 = 100
    nmeas::Int64 = 1
    thermalized::Bool = false
end
export HMC

struct FallbackHMC <: AbstractHMC
    params::HMCParams
end

sampler(lftws::AbstractLFT, hmcp::HMCParams) = FallbackHMC(hmcp)
sample!(lftws::AbstractLFT, samplerws::AbstractHMC, debugger::Union{Vector{<:AbstractDebugger},Nothing} = nothing) = hmc!(lftws, samplerws, debugger)
