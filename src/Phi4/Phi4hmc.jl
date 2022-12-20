
function generate_momenta!(phiws::Phi4, lp::Phi4Parm)
    # Create momenta for phi
    phiws.mom .= Random.randn(size(phiws.phi))
end

function Hamiltonian(phiws::Phi4, lp::Phi4Parm)
    return Hamiltonian(phiws.mom, phiws, lp)
end

function Hamiltonian(mom, phiws::Phi4, lp::Phi4Parm)
    H = mapreduce(x -> x^2, +, mom)/2.0 + action(phiws, lp)
    return H
end

function update_momenta!(phiws::Phi4, epsilon, lp::Phi4Parm)
    update_momenta!(phiws.mom, phiws, epsilon, lp)
    return nothing
end

function update_momenta!(mom, phiws::Phi4, epsilon, lp::Phi4Parm)

    # Load phi force
    force!(phiws, lp) 

    # Update phi momenta
    update_momenta!(mom, phiws.frc, epsilon, phiws, lp)

    return nothing
end

function update_momenta!(mom, frc, epsilon, phiws::Phi4, lp::Phi4Parm)
    mom .= mom .+ epsilon .* frc
    return nothing
end

function update_fields!(phiws::Phi4, epsilon, lp::Phi4Parm)
    update_fields!(phiws.phi, phiws.mom, epsilon, phiws, lp)
    return nothing
end

function update_fields!(phi, mom, epsilon, phiws::Phi4, lp::Phi4Parm) 
    # Update phi field
    phi .= phi .+ epsilon * mom
    return nothing
end

