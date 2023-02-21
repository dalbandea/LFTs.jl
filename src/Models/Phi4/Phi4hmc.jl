
function generate_momenta!(phiws::Phi4, hmcws::Phi4HMC)
    # Create momenta for phi
    hmcws.mom .= Random.randn(size(phiws.phi))
end

function Hamiltonian(phiws::Phi4, hmcws::Phi4HMC)
    H = mapreduce(x -> x^2, +, hmcws.mom)/2.0 + action(phiws)
    return H
end

function update_momenta!(phiws::Phi4, epsilon, hmcws::Phi4HMC)

    # Load phi force
    force!(phiws, hmcws) 

    # Update phi momenta
    hmcws.mom .= hmcws.mom .+ epsilon .* hmcws.frc

    return nothing
end

function update_fields!(phiws::Phi4, epsilon, hmcws::Phi4HMC)
    # Update phi field
    phiws.phi .= phiws.phi .+ epsilon .* hmcws.mom
    return nothing
end
