
function action(phiws::Phi4workspace, lp::Phi4Parm) where {T}
    return action(phiws.phi, phiws, lp)
end

function action(phi, phiws::Phi4workspace, lp::Phi4Parm) where {T}

    # For type stability initialize to correct type
    # https://www.juliabloggers.com/writing-type-stable-julia-code/
    act = zero(eltype(phi))

    Nx = lp.iL[1]
    Ny = lp.iL[2]
    for i in 1:Nx, j in 1:Ny
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        act = act - lp.beta * phi[i,j] * (phi[iu,j] + phi[i,ju]) +
                (1 - 2 * lp.lambda) * phi[i,j]^2 + lp.lambda * phi[i,j]^4
    end

    return act
end


function force!(phiws::Phi4workspace, lp::Phi4Parm) 
    return force!(phiws.frc, phiws.phi, phiws, lp)
end

# Force is computed in place for performance
function force!(frc, phi, phiws::Phi4workspace, lp::Phi4Parm) 

    Nx = size(phi,1)
    Ny = size(phi,2)
    for i in 1:Nx, j in 1:Ny
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        id = mod1(i-1, Nx)
        jd = mod1(j-1, Ny)
        
        frc[i,j] = lp.beta * (phi[iu,j] + phi[id,j] + phi[i,ju] + phi[i,jd]) + 
                    2 * phi[i,j] * (2 * lp.lambda * (1 - phi[i,j]^2) - 1)
    end

    return nothing
end

