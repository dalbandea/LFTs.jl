
function check_force(phiws::Phi4, hmcws::Phi4HMC, epsilon)
    phiws2 = deepcopy(phiws)

    # Initial action
    Si = action(phiws)

    Nx = phiws.params.iL[1]
    Ny = phiws.params.iL[2]
    V = Nx * Ny

    # Compute force of configuration
    force!(phiws, hmcws)

    F_diff = 0.0

    for i in 1:Nx, j in 1:Ny
        phiws2.phi[j,i] += epsilon

        # Final action
        Sf = action(phiws2)

        # Numerical force
        F_num = (Sf - Si)/epsilon

        # Analytical force
        F_ana = hmcws.frc[j,i]

        # Difference
        F_diff = F_ana + F_num

        phiws2.phi[j,i] -= epsilon
    end

    println("Difference: ", F_diff / V)
end


function reversibility!(phiws::Phi4, hmcws::Phi4HMC, epsilon, ns)
    phi_cp = copy(phiws.phi)

    # Create momenta
    generate_momenta!(phiws, hmcws)

    leapfrog!(phiws, hmcws, epsilon, ns)
    hmcws.mom .= -hmcws.mom
    leapfrog!(phiws, hmcws, epsilon, ns)

    dphi2 = sum(abs2, phiws.phi .- phi_cp)

    println("|Δϕ|² = ", dphi2)

    return nothing
end
