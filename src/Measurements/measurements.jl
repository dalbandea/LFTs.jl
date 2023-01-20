abstract type AbstractObservable end
abstract type AbstractScalar <: AbstractObservable end
abstract type AbstractCorrelator <: AbstractObservable end

function measure(observable::AbstractObservable, lftws::LFTworkspace, lp::LattParm) end
function write(observable::AbstractObservable) end
function save(observable::AbstractObservable) end
function read(observable::AbstractObservable) end
function analyze(observable::AbstractObservable) end

function measure(observables::Array{T}, lftws::LFTworkspace, lp::LattParm) where T <: AbstractObservable 
    for observable in observables
        measure(observable, lftws, lp)
    end
end 

function measure(observable::AbstractScalar, lftws::LFTworkspace, lp::LattParm)
    observable.result = observable(lftws, lp)
    return nothing
end
export measure

function write(obs::AbstractScalar)
    global io_stat = open(obs.filepath, "a")
    write(io_stat, "$(obs.result)\n")
    close(io_stat)
    return nothing
end

save(obs::AbstractScalar) = push!(obs.history, obs.result)

read(obs::AbstractScalar) = vec(DelimitedFiles.readdlm(obs.filepath))



function sample_and_measure!(observables::Array{T}, lftws::LFTworkspace, sampler::AbstractSampler, lp::LattParm; verbose::Bool = false, get_history::Bool = false) where T <: AbstractObservable

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
            for observable in observables
                measure(observable, lftws, lp)
                get_history ? save(observable) : write(observable)
            end
        end
    end
end
export sample_and_measure!

