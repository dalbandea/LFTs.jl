abstract type Integrator end

struct Leapfrog <: Integrator 
    epsilon::Float64
    nsteps::Int64
    function Leapfrog(tau::Float64, nsteps::Int64)
        epsilon = tau / nsteps
        return new(epsilon, nsteps)
    end
end

struct OMF4 <: Integrator
    epsilon::Float64
    nsteps::Int64
    function OMF4(tau::Float64, nsteps::Int64)
        epsilon = tau / nsteps
        return new(epsilon, nsteps)
    end
end



