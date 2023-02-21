
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


struct QuantumRotorWorkspace{T, N, P <: LFTParm} <: QuantumRotor
    PRC::Type{T}
    phi::Array{T, N}
    params::P
    function QuantumRotorWorkspace(::Type{T}, lp::QuantumRotorParm) where {T <: AbstractFloat}
        phi = Array{T, 1}(undef, lp.iT)
        return new{T, 1, typeof(lp)}(T, phi, lp)
    end
end

function (::Type{QuantumRotor})(::Type{T} = Float64; kwargs...) where {T <: AbstractFloat}
    return QuantumRotorWorkspace(T, QuantumRotorParm(;kwargs...))
end



struct QuantumRotorHMC{A <: AbstractArray} <: AbstractHMC
    params::HMCParams
    frc::A
    mom::A
end

function QuantumRotorHMC(phiws::QuantumRotor, hmcp::HMCParams)
    frc = similar(phiws.phi)
    mom = similar(phiws.phi)
    return QuantumRotorHMC{typeof(frc)}(hmcp, frc, mom)
end

sampler(lftws::QuantumRotor, hmcp::HMCParams) = QuantumRotorHMC(lftws, hmcp)

function copy!(phiws_dst::QuantumRotor, phiws_src::QuantumRotor, hmcws::QuantumRotorHMC)
    phiws_dst.phi .= phiws_src.phi
    return nothing
end

function randomize!(phiws::QuantumRotor)
    phiws.phi .= 2pi*Random.rand(phiws.PRC, size(phiws.phi)...) .- pi
    return nothing
end



#custom mod
#changes mod domain from (0 to z) to (-z/2 to z/2)
function Mod(x, z)

    a = mod(x,z)
    if a <= z/2
        return a
    else
        return a - z
    end

end

function action(qrws::QuantumRotor)
    return action(qrws, qrws.params.disc, qrws.params.BC)
end

function action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{B}) where B <: AbstractBoundaryCondition
    S = zero(qrws.PRC)

    for t in 1:qrws.params.iT-1
        S += Mod(qrws.phi[t+1] - qrws.phi[t], 2pi)^2
    end

    S += boundary_action(qrws, disc, BC)

    return qrws.params.I * S / 2
end

boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{PeriodicBC}) = Mod(qrws.phi[1] - qrws.phi[end], 2pi)^2
boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{AntiperiodicBC}) = Mod(-qrws.phi[1] - qrws.phi[end], 2pi)^2
boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{OpenBC}) = zero(qrws.PRC)


function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    return force!(qrws, hmcws, qrws.params.disc, qrws.params.BC)
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{ClassicalPerfectDiscretization}, BC::Type{B}) where B <: AbstractBoundaryCondition

    for t in 2:qrws.params.iT-1
        hmcws.frc[t] = -qrws.params.I * (Mod(qrws.phi[t]-qrws.phi[t-1], 2pi) - Mod(qrws.phi[t+1] - qrws.phi[t], 2pi))
    end

    boundary_force!(qrws, hmcws, disc, BC)
    
    return nothing
end


function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{ClassicalPerfectDiscretization}, BC::Type{PeriodicBC})
    hmcws.frc[1] = -qrws.params.I * (Mod(qrws.phi[1] - qrws.phi[end], 2pi) - Mod(qrws.phi[2] - qrws.phi[1], 2pi))
    hmcws.frc[end] = -qrws.params.I * (Mod(qrws.phi[end] - qrws.phi[end-1], 2pi) - Mod(qrws.phi[1] - qrws.phi[end], 2pi))
    return nothing
end

function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    for i in 1:length(hmcws.mom)
        hmcws.mom[i] = randn()
    end
    return nothing
end

function Hamiltonian(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    H = mapreduce(x -> x^2, +, hmcws.mom) + action(qrws)
    return H
end

function update_momenta!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC)

    # Load phi force
    force!(qrws, hmcws) 

    # Update phi momenta
    hmcws.mom .= hmcws.mom .+ epsilon .* hmcws.frc

    return nothing
end


function update_fields!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC)
    # Update phi field
    qrws.phi .= qrws.phi .+ epsilon .* hmcws.mom
    return nothing
end

function flip_momenta_sign!(hmcws::QuantumRotorHMC)
    hmcws.mom .= .- hmcws.mom
    return nothing
end
