abstract type AbstractSampler end
abstract type AbstractHMC <: AbstractSampler end
abstract type AbstractObservable end

function sample!(lftws::LFTworkspace, sampler::AbstractSampler, lp::LattParm) end 
export sample!
function measure(observable::AbstractObservable, lftws::LFTworkspace, lp::LattParm) end
function measure(observables::Array{T}, lftws::LFTworkspace, lp::LattParm) where T <: AbstractObservable 
    for observable in observables
        measure(observable, lftws, lp)
    end
end 
function analyze(observable::AbstractObservable) end
export analyze


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

function sample_and_measure!(observables::Array{T}, lftws::LFTworkspace, sampler::AbstractSampler, lp::LattParm; verbose::Bool = false) where T <: AbstractObservable

    if sampler.thermalized == false
        for i in 1:sampler.ntherm
            verbose && print("Thermalizing... $(i)\r") 
            sample!(lftws, sampler, lp)
        end
        sampler.thermalized == true
    end

    for i in 1:sampler.ntraj
        verbose && print("Sampling... $(i)\r") 
        sample!(lftws, sampler, lp)

        if i%sampler.nmeas == 0
            verbose && print("Measuring...\r")
            measure(observables, lftws, lp)
        end
    end
end
export sample_and_measure!

