function sample!(qrws::QuantumRotor, mws::AbstractMetropolis; do_winding = false)
    sweep!(qrws, mws)
    do_winding && winding_step!(qrws)
    return nothing
end

function sweep!(qrws::QuantumRotor, mws::AbstractMetropolis)
    for t in 1:qrws.params.iT
        update_variable!(qrws, mws, t)
    end
    return nothing
end

function update_variable!(qrws::QuantumRotor, mws::AbstractMetropolis, t::Int64)

    ws_cp = deepcopy(qrws)

    #change in phi
    delta = mws.params.weight * randn()

    sini = action(qrws)

    qrws.phi[t] += delta

    sfin = action(qrws)

    ds = sfin - sini

    metropolis_accept_reject!(qrws, ws_cp, ds)

    return nothing
end

function winding_step!(qrws::QuantumRotor)
    ws_cp = deepcopy(qrws)

    sini = action(qrws)

    r = rand()

    if r > 0.5
        winding!(qrws)
    else
        antiwinding!(qrws)
    end

    sfin = action(qrws)

    ds = sfin - sini

    metropolis_accept_reject!(qrws, ws_cp, ds)

    return nothing
end
