
function action_check(phiws::Phi4Upscaledworkspace, lp::Phi4Parm)
    phi_interp = interpolation_upscale(phiws.phi)
    Nx = size(phi_interp, 1)
    Ny = size(phi_interp, 2)

    S_normal = phi4action(phi_interp, Nx, Ny, lp.beta, lp.lambda)
    S_integrated_out = action(phiws, lp)

    println("Action without integrating out:    ", S_normal)
    println("Action after integrating out:      ", S_integrated_out)
end

function interpolation_check(phiws::Phi4Upscaledworkspace)
    phi_interp = interpolation_upscale(phiws.phi)

    # Check that sum over upscaled config is 4 times the sum over the original
    println("Diff: ", sum(phi_interp) - 4*sum(phiws.phi))
    return nothing
end
