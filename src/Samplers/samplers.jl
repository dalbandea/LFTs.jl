abstract type AbstractSampler end
abstract type AbstractHMC <: AbstractSampler end

function sample!(lftws::LFTworkspace, sampler::AbstractSampler, lp::LattParm) end 
export sample!

include("HMC/integrators/integrators.jl")
export Leapfrog, OMF4
include("HMC/hmc.jl")
export hmc!

Base.@kwdef mutable struct HMC <: AbstractHMC
    integrator::Integrator = Leapfrog()
    ntherm::Int64 = 10
    ntraj::Int64 = 100
    nmeas::Int64 = 1
    thermalized::Bool = false
end
export HMC

sample!(lftws::LFTworkspace, sampler::HMC, lp::LattParm) = hmc!(lftws, sampler.integrator, lp)

