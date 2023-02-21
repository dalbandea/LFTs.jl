abstract type AbstractObservable end
abstract type AbstractScalar <: AbstractObservable end
abstract type AbstractCorrelator <: AbstractObservable end

function initialize_sampler(sampler::AbstractSampler, lftws::AbstractLFT) end

measure!(observable::AbstractObservable, lftws::AbstractLFT) = observable(lftws)
function write(observable::AbstractObservable) end
function save!(observable::AbstractObservable) end
function read(observable::AbstractObservable) end
function analyze(observable::AbstractObservable) end

function measure!(observables::Array{T}, lftws::AbstractLFT) where T <: AbstractObservable 
    for observable in observables
        measure!(observable, lftws)
    end
end 

export measure!

## AbstractScalar

function write(obs::AbstractScalar)
    global io_stat = open(obs.filepath, "a")
    write(io_stat, "$(obs.result)\n")
    close(io_stat)
    return nothing
end

function save!(obs::AbstractScalar)
    push!(obs.history, obs.result)
    return nothing
end

read(obs::AbstractScalar) = vec(DelimitedFiles.readdlm(obs.filepath))

function read(::Type{T}, filepath::String; column::Int64 = 1) where T <: AbstractScalar
    return DelimitedFiles.readdlm(filepath, ',', Float64)[:,1]
end

## AbstractCorrelator

function write(obs::AbstractCorrelator)
    global io_stat = open(obs.filepath, "a")
    write(io_stat, "$(obs.result[1])")
    for i in 2:length(obs.result)
        write(io_stat, ",$(obs.result[i])")
    end
    write(io_stat, "\n")
    close(io_stat)
    return nothing
end

function save!(obs::AbstractCorrelator)
    push!(obs.history, obs.result)
    return nothing
end
export save!

function read(::Type{T}, filepath::String) where T <: AbstractCorrelator
    return DelimitedFiles.readdlm(filepath, ',', Float64)
end


## General sampling and measure method

function sample_and_measure!(observables::Array{T}, lftws::AbstractLFT, sampler::AbstractSampler; verbose::Bool = false, get_history::Bool = false) where T <: AbstractObservable

    if sampler.thermalized == false
        for i in 1:sampler.ntherm
            verbose && print("Thermalizing... $(i)\r") 
            sample!(lftws, sampler)
        end
        sampler.thermalized == true
    end

    for i in 1:sampler.ntraj
        verbose && print("Sampling... $(i)\r") 
        sample!(lftws, sampler)

        if i%sampler.nmeas == 0
            verbose && print("Measuring...\r")
            for observable in observables
                measure!(observable, lftws)
                get_history ? save!(observable) : write(observable)
            end
        end
    end
end
export sample_and_measure!

