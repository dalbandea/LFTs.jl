function action(phiws::Phi4Upscaledworkspace, lp::Phi4Parm)
    return action(phiws.phi, phiws, lp)
end

function action(phi, phiws::Phi4Upscaledworkspace, lp::Phi4Parm) where {T}

    # For type stability initialize to correct type
    # https://www.juliabloggers.com/writing-type-stable-julia-code/
    act = zero(eltype(phi))

    Nx = lp.iL[1]
    Ny = lp.iL[2]
    for j in 1:Ny, i in 1:Nx
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        act = act - 
        lp.beta * ( 1/2 * (phi[i,j] + phi[iu,j])^2 + 1/2*(phi[i,j] + phi[i,ju])^2 +
                   1/4 * plaquette(phiws, i, j, lp)^2 ) + 
        (1 - 2 * lp.lambda) * (phi[i,j]^2 + 1/4*(phi[i,j] + phi[iu,j])^2 +
                               1/4*(phi[i,j] + phi[i,ju])^2 + 1/16 *
                               plaquette(phiws,i,j,lp)^2) + 
        lp.lambda * (phi[i,j]^4 + 1/16*(phi[i,j] + phi[iu,j])^4 + 1/16*(phi[i,j] 
                        + phi[i,ju])^4 + 1/256 * plaquette(phiws, i, j,lp)^4) 
    end

    return act
end

function plaquette(phiws::Phi4Upscaledworkspace, i::Int64, j::Int64, lp::Phi4Parm)
    return plaquette(phiws.phi, i, j, lp.iL[1], lp.iL[2], phiws)
end


function plaquette(phi, i, j, Nx, Ny, phiws::Phi4Upscaledworkspace)
    iu = mod1(i+1, Nx)
    ju = mod1(j+1, Ny)

    return phi[i,j]+phi[iu,j]+phi[i,ju]+phi[iu,ju]
end

function force!(phiws::Phi4Upscaledworkspace, lp::Phi4Parm) 
    return force!(phiws.frc, phiws.phi, phiws, lp)
end

# Force is computed in place for performance
function force!(frc, phi, phiws::Phi4Upscaledworkspace, lp::Phi4Parm) 

    Nx = size(phi,1)
    Ny = size(phi,2)
    for j in 1:Ny, i in 1:Nx
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        id = mod1(i-1, Nx)
        jd = mod1(j-1, Ny)

        frc[i,j] = lp.beta *  
        ( (4 * phi[i,j] + phi[iu,j] + phi[id,j] + phi[i,ju] + phi[i,jd])  + 
         1/2 * (plaquette(phiws,i,j,lp)+plaquette(phiws,id,j,lp)+plaquette(phiws,i,jd,lp)+plaquette(phiws,id,jd,lp)) )  -
        (1 - 2 * lp.lambda) * 
        (2 * phi[i,j] + 
         1/2*(4 * phi[i,j] + phi[iu,j] + phi[id,j] + phi[i,ju] + phi[i,jd]) +
         1/8*(plaquette(phiws,i,j,lp)+plaquette(phiws,id,j,lp)+plaquette(phiws,i,jd,lp)+plaquette(phiws,id,jd,lp))) -
        lp.lambda * 
        (4 * phi[i,j]^3 + 
         1/4 * ((phi[i,j] + phi[iu,j])^3 + (phi[i,j] + phi[id,j])^3 + (phi[i,j] + phi[i,ju])^3 + (phi[i,j] + phi[i,jd])^3 ) +
         1/64 * ( plaquette(phiws,i,j,lp)^3 + plaquette(phiws,id,j,lp)^3 + plaquette(phiws,i,jd,lp)^3 + plaquette(phiws,id,jd,lp)^3 )) 
    end

    return nothing
end

