
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



