abstract type AbstractIntegrator end

Base.@kwdef struct Leapfrog <: AbstractIntegrator 
    epsilon::Float64 = 0.1 
    nsteps::Int64 = 10
    function Leapfrog(tau::Float64, nsteps::Int64)
        epsilon = tau / nsteps
        return new(epsilon, nsteps)
    end
end

struct OMF4 <: AbstractIntegrator
    epsilon::Float64
    nsteps::Int64
    function OMF4(tau::Float64, nsteps::Int64)
        epsilon = tau / nsteps
        return new(epsilon, nsteps)
    end
end



