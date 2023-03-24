abstract type AbstractSampler end
abstract type SamplerParameters end

abstract type AbstractHMC <: AbstractSampler end
abstract type HMCParams <: SamplerParameters end

abstract type AbstractHMCLFT end


abstract type AbstractMetropolis <: AbstractSampler end
abstract type MetropolisParams <: SamplerParameters end


function copy!(lftws_dest::L, lftws_src::L) where L <: AbstractLFT
    error("No function copy! for $(typeof(lftws_dest))")
    return nothing
end

function sampler(lftws::AbstractLFT, samplerparms::SamplerParameters) end
function sample!(lftws::AbstractLFT, samplerws::AbstractSampler) end 
function sample!(lftws::AbstractLFT, samplerparms::SamplerParameters)
    samplerws = sampler(lftws, samplerparms)
    return sample!(lftws, samplerws)
end
export sample!

function metropolis_accept_reject!(lftws::L, lftcp::L, samplerws::S, dS::Float64) where {L <: AbstractLFT, S <: AbstractSampler}
    pacc = exp(-dS)
    if (pacc < 1.0)
        r = rand()
        if (r > pacc) 
            copy!(lftws, lftcp)
            @info("    REJECT: Energy [difference: $(dS)]")
        else
            @info("    ACCEPT:  Energy [difference: $(dS)]")
        end
    else
        @info("    ACCEPT:  Energy [difference: $(dS)]")
    end
    return nothing
end


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

include("Metropolis/metropolis.jl")

Base.@kwdef mutable struct Metropolis <: MetropolisParams
    weight::Float64 = 0.1
    ntherm::Int64 = 10
    ntraj::Int64 = 100
    nmeas::Int64 = 1
    thermalized::Bool = false
end
export Metropolis

struct FallbackMetropolis{P <: MetropolisParams} <: AbstractMetropolis
    params::P
end

sampler(lftws::AbstractLFT, hmcp::HMCParams) = FallbackHMC(hmcp)
sample!(lftws::AbstractLFT, samplerws::AbstractHMC) = hmc!(lftws, samplerws)

sampler(lftws::AbstractLFT, mp::MetropolisParams) = FallbackMetropolis(mp)
sample!(lftws::AbstractLFT, samplerws::AbstractMetropolis) = sweep!(lftws, samplerws)
