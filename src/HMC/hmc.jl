
# These functions need to be defined for every model (subtype of LFTworkspace)
function generate_momenta!(lftws::LFTworkspace, lp::LattParm) end
function Hamiltonian(lftws::LFTworkspace, lp::LattParm) end
function copy!(lftws_dest::LFTworkspace, lftws_src::LFTworkspace, lp::LattParm) end
function update_momenta!(lftws::LFTworkspace, epsilon, lp::LattParm) end
function update_fields!(lftws::LFTworkspace, epsilon, lp::LattParm) end

function HMC!(lftws::LFTworkspace, epsilon, ns, lp::LattParm)
    # Create copy of current configuration
    ws_cp = deepcopy(lftws)

    # Generate random momenta
    generate_momenta!(lftws, lp)

    # Compute initial Hamiltonian
    hini = Hamiltonian(lftws, lp)

    # Molecular Dynamics
    leapfrog!(lftws, epsilon, ns, lp)

    # Compute final Hamiltonian
    hfin = Hamiltonian(lftws, lp)

    dH = hfin - hini
    pacc = exp(-dH)
    if (pacc < 1.0)
        r = rand()
        if (r > pacc) 
            copy!(lftws, ws_cp, lp)
            @info("    REJECT: Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
        else
            @info("    ACCEPT:  Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
        end
    else
        @info("    ACCEPT:  Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
    end

    return dH
end

function leapfrog!(lftws::LFTworkspace, epsilon, ns, lp::LattParm)

	# First half-step for momenta
    update_momenta!(lftws, epsilon/2.0, lp)

	# ns-1 steps
	for i in 1:(ns-1) 
		# Update fields
        update_fields!(lftws, epsilon, lp) 

		#Update momenta
        update_momenta!(lftws, epsilon, lp)
	end
	# Last update for fields
    update_fields!(lftws, epsilon, lp) 

	# Last half-step for momenta
    update_momenta!(lftws, epsilon/2.0, lp)

	return nothing
end

