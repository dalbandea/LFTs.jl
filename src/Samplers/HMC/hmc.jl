
# These functions need to be defined for every model (subtype of LFTworkspace)
## Mandatory
function generate_momenta!(lftws::LFTworkspace, lp::LattParm) end
function Hamiltonian(lftws::LFTworkspace, lp::LattParm) end
function action(lftws::LFTworkspace, lp::LattParm) end
function copy!(lftws_dest::LFTworkspace, lftws_src::LFTworkspace, lp::LattParm) end
function update_momenta!(lftws::LFTworkspace, epsilon, lp::LattParm) end
function update_fields!(lftws::LFTworkspace, epsilon, lp::LattParm) end

## Optional
function generate_pseudofermions!(lftws::LFTworkspace, lp::LattParm) end


molecular_dynamics!(lftws::LFTworkspace, integr::Leapfrog, lp::LattParm) =
                        leapfrog!(lftws, integr.epsilon, integr.nsteps, lp)
molecular_dynamics!(lftws::LFTworkspace, integr::OMF4, lp::LattParm) =
                        OMF4!(lftws, integr.epsilon, integr.nsteps, lp)

function HMC!(lftws::LFTworkspace, integrator::Integrator, lp::LattParm)
    # Create copy of current configuration
    ws_cp = deepcopy(lftws)

    # Generate random momenta
    generate_momenta!(lftws, lp)

    # Initialize pseudofermion and related fields
    generate_pseudofermions!(lftws, lp)

    # Compute initial Hamiltonian
    hini = Hamiltonian(lftws, lp)

    # Molecular Dynamics
    molecular_dynamics!(lftws, integrator, lp)

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



function OMF4!(lftws::LFTworkspace, epsilon, ns, lp::LattParm)

    r1::Float64 =  0.08398315262876693
    r2::Float64 =  0.2539785108410595
    r3::Float64 =  0.6822365335719091
    r4::Float64 = -0.03230286765269967
    r5::Float64 =  0.5-r1-r3
    r6::Float64 =  1.0-2.0*(r2+r4)

    frc = zeros(Float64, lp.iL[1], lp.iL[2])

    for i in 1:ns
        # STEP 1
        update_momenta!(lftws, r1*epsilon, lp)
        update_fields!(lftws, r2*epsilon, lp) 

        # STEP 2
        update_momenta!(lftws, r3*epsilon, lp)
        update_fields!(lftws, r4*epsilon, lp) 

        # STEP 3
        update_momenta!(lftws, r5*epsilon, lp)
        update_fields!(lftws, r6*epsilon, lp) 

        # STEP 4
        update_momenta!(lftws, r5*epsilon, lp)
        update_fields!(lftws, r4*epsilon, lp) 

        # STEP 5
        update_momenta!(lftws, r3*epsilon, lp)
        update_fields!(lftws, r2*epsilon, lp) 

        # STEP 6
        update_momenta!(lftws, r1*epsilon, lp)
    end

    return nothing
end

