
function action(phiws::Phi4, lp::Phi4Parm) where {T}
    return phi4action(phiws.phi, lp.iL[1], lp.iL[2], lp.beta, lp.lambda)
end

function phi4action(phi, Nx, Ny, beta, lambda)
    # For type stability initialize to correct type
    # https://www.juliabloggers.com/writing-type-stable-julia-code/
    act = zero(eltype(phi))

    for j in 1:Ny, i in 1:Nx
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        act = act - beta * phi[i,j] * (phi[iu,j] + phi[i,ju]) +
                (1 - 2 * lambda) * phi[i,j]^2 + lambda * phi[i,j]^4
    end

    return act
end

function phi4action(phi, beta, lambda)
    Nx = size(phi, 1)
    Ny = size(phi, 2)
    return phi4action(phi, Nx, Ny, beta, lambda)
end

function force!(phiws::Phi4, lp::Phi4Parm) 
    return phi4force!(phiws.frc, phiws.phi, lp.iL[1], lp.iL[2], lp.beta, lp.lambda)
end

# Force is computed in place for performance
function phi4force!(frc, phi, Nx, Ny, beta, lambda) 

    for j in 1:Ny, i in 1:Nx
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        id = mod1(i-1, Nx)
        jd = mod1(j-1, Ny)
        
        frc[i,j] = beta * (phi[iu,j] + phi[id,j] + phi[i,ju] + phi[i,jd]) + 
                    2 * phi[i,j] * (2 * lambda * (1 - phi[i,j]^2) - 1)
    end

    return nothing
end

