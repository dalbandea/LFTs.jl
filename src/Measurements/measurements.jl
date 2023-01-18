abstract type AbstractObservable end

function measure(observable::AbstractObservable, lftws::LFTworkspace, lp::LattParm) end
function measure(observables::Array{T}, lftws::LFTworkspace, lp::LattParm) where T <: AbstractObservable 
    for observable in observables
        measure(observable, lftws, lp)
    end
end 
function analyze(observable::AbstractObservable) end
export measure, analyze

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

