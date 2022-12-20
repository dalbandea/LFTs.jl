
function check_force(phiws::Phi4, lp::Phi4Parm, epsilon)
    phiws2 = deepcopy(phiws)

    # Initial action
    Si = action(phiws, lp)

    V = lp.iL[1] * lp.iL[2]
    Nx = lp.iL[1]
    Ny = lp.iL[2]

    # Compute force of configuration
    force!(phiws, lp)

    F_diff = 0.0

    for i in 1:Nx, j in 1:Ny
        phiws2.phi[j,i] += epsilon

        # Final action
        Sf = action(phiws2, lp)

        # Numerical force
        F_num = (Sf - Si)/epsilon

        # Analytical force
        F_ana = phiws.frc[j,i]

        # Difference
        F_diff = F_ana + F_num

        phiws2.phi[j,i] -= epsilon
    end

    println("Difference: ", F_diff / V)
end


function reversibility!(phiws::Phi4, epsilon, ns, lp::LattParm)
    phi_cp = copy(phiws.phi)

    # Create momenta
    generate_momenta!(phiws, lp)

    leapfrog!(phiws, epsilon, ns, lp)
    phiws.mom .= -phiws.mom
    leapfrog!(phiws, epsilon, ns, lp)

    dphi2 = sum(abs2, phiws.phi .- phi_cp)

    println("|Δϕ|² = ", dphi2)

    return nothing
end
