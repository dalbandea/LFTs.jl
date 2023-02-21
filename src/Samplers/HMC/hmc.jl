# These functions need to be defined for every model (subtype of AbstractLFT)
## Mandatory
function generate_momenta!(lftws::AbstractLFT, hmcws::AbstractHMC) 
    error("No function generate_momenta! for $(typeof(lftws))")
    return nothing
end
function Hamiltonian(lftws::AbstractLFT, hmcws::AbstractHMC) 
    error("No function Hamiltonian for $(typeof(lftws))")
    return nothing
end
function action(lftws::AbstractLFT) 
    error("No function action for $(typeof(lftws))")
    return nothing
end
action(lftws::AbstractLFT, hmcws::AbstractHMC) = action(lftws)
function copy!(lftws_dest::AbstractLFT, lftws_src::AbstractLFT, hmcws::AbstractHMC)
    error("No function copy! for $(typeof(lftws_dest))")
    return nothing
end
function update_momenta!(lftws::AbstractLFT, epsilon, hmcws::AbstractHMC) 
    error("No function update_momenta! for $(typeof(lftws))")
    return nothing
end
function update_fields!(lftws::AbstractLFT, epsilon, hmcws::AbstractHMC) 
    error("No function update_fields! for $(typeof(lftws))")
    return nothing
end

## Optional
function generate_pseudofermions!(lftws::AbstractLFT, hmcws::AbstractHMC) end

molecular_dynamics!(lftws::AbstractLFT, hmcws::AbstractHMC) =
                        molecular_dynamics!(lftws, hmcws, hmcws.params.integrator)
molecular_dynamics!(lftws::AbstractLFT, hmcws::AbstractHMC, integr::Leapfrog) =
                        leapfrog!(lftws, hmcws, integr.epsilon, integr.nsteps)
molecular_dynamics!(lftws::AbstractLFT, hmcws::AbstractHMC, integr::OMF4) =
                        OMF4!(lftws, hmcws, integr.epsilon, integr.nsteps)

function hmc!(lftws::AbstractLFT, hmcws::AbstractHMC)
    # Create copy of current configuration
    ws_cp = deepcopy(lftws)

    # Generate random momenta
    generate_momenta!(lftws, hmcws)

    # Initialize pseudofermion and related fields
    generate_pseudofermions!(lftws, hmcws)

    # Compute initial Hamiltonian
    hini = Hamiltonian(lftws, hmcws)

    # Molecular Dynamics
    molecular_dynamics!(lftws, hmcws)

    # Compute final Hamiltonian
    hfin = Hamiltonian(lftws, hmcws)

    dH = hfin - hini
    pacc = exp(-dH)
    if (pacc < 1.0)
        r = rand()
        if (r > pacc) 
            copy!(lftws, ws_cp, hmcws)
            @info("    REJECT: Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
        else
            @info("    ACCEPT:  Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
        end
    else
        @info("    ACCEPT:  Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]")
    end

    return dH
end

function leapfrog!(lftws::AbstractLFT, hmcws::AbstractHMC, epsilon, ns)

	# First half-step for momenta
    update_momenta!(lftws, epsilon/2.0, hmcws)

	# ns-1 steps
	for i in 1:(ns-1) 
		# Update fields
        update_fields!(lftws, epsilon, hmcws) 

		#Update momenta
        update_momenta!(lftws, epsilon, hmcws)
	end
	# Last update for fields
    update_fields!(lftws, epsilon, hmcws) 

	# Last half-step for momenta
    update_momenta!(lftws, epsilon/2.0, hmcws)

	return nothing
end



function OMF4!(lftws::AbstractLFT, hmcws::AbstractHMC, epsilon, ns)

    r1::Float64 =  0.08398315262876693
    r2::Float64 =  0.2539785108410595
    r3::Float64 =  0.6822365335719091
    r4::Float64 = -0.03230286765269967
    r5::Float64 =  0.5-r1-r3
    r6::Float64 =  1.0-2.0*(r2+r4)

    frc = zeros(Float64, lftws.iL[1], lftws.iL[2])

    for i in 1:ns
        # STEP 1
        update_momenta!(lftws, r1*epsilon, hmcws)
        update_fields!(lftws, r2*epsilon, hmcws) 

        # STEP 2
        update_momenta!(lftws, r3*epsilon, hmcws)
        update_fields!(lftws, r4*epsilon, hmcws) 

        # STEP 3
        update_momenta!(lftws, r5*epsilon, hmcws)
        update_fields!(lftws, r6*epsilon, hmcws) 

        # STEP 4
        update_momenta!(lftws, r5*epsilon, hmcws)
        update_fields!(lftws, r4*epsilon, hmcws) 

        # STEP 5
        update_momenta!(lftws, r3*epsilon, hmcws)
        update_fields!(lftws, r2*epsilon, hmcws) 

        # STEP 6
        update_momenta!(lftws, r1*epsilon, hmcws)
    end

    return nothing
end

