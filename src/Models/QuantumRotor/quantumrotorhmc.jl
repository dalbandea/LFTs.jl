
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
