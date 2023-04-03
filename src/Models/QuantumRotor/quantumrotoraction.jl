
function action(qrws::QuantumRotor)
    S = zero(qrws.PRC)

    for t in 1:qrws.params.iT
        S += action_t(qrws, t)
    end

    return qrws.params.I * S / 2
end

action_t(qrws::QuantumRotor, t::Int64) = action_t(qrws, qrws.params.disc, t)
function action_t(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, t::Int64)
    tu = right(qrws, t)
    ds = Mod(qrws.phi[tu] - qrws.phi[t], 2pi)^2
    return ds
end

"""
    daction_t

Compute action of points coupled to the angle `t`.
"""
daction_t(qrws::QuantumRotor, t::Int64) = daction_t(qrws, qrws.params.disc, t)
function daction_t(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, t::Int64)
    td = left(qrws, t)
    ds = action_t(qrws, td) + action_t(qrws, t)
    return ds
end

right(qrws::QuantumRotor, t) = right(qrws, qrws.params.BC, t)
function right(qrws::QuantumRotor, BC::Type{B}, t::Int64) where B <: AbstractBoundaryCondition 
    error("Function right not implemented for boundary condition $B")
end
right(qrws::QuantumRotor, BC::Type{PeriodicBC}, t::Int64) = mod1(t+1,
                                                                 qrws.params.iT)

left(qrws::QuantumRotor, t) = left(qrws, qrws.params.BC, t)
function left(qrws::QuantumRotor, BC::Type{B}, t::Int64) where B <: AbstractBoundaryCondition 
    error("Function left not implemented for boundary condition $B")
end
left(qrws::QuantumRotor, BC::Type{PeriodicBC}, t::Int64) = mod1(t-1,
                                                                 qrws.params.iT)

boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{PeriodicBC}) = Mod(qrws.phi[1] - qrws.phi[end], 2pi)^2
boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{AntiperiodicBC}) = Mod(-qrws.phi[1] - qrws.phi[end], 2pi)^2
boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{OpenBC}) = zero(qrws.PRC)



