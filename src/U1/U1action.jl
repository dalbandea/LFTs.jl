# function U1plaquette!(plx, U, Nx, Ny)

#     i1 = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#     i2 = (CUDA.blockIdx().y - 1) * CUDA.blockDim().y + CUDA.threadIdx().y

#     iu1 = mod(i1, Nx) + 1
#     iu2 = mod(i2, Ny) + 1
    
#     plx[i1,i2] = real(U[i1,i2,1] *
#                       U[iu1,i2,2] *
#                       conj(U[i1,iu2,1] *
#                            U[i1,i2,2]))
    
#     return nothing
# end

KernelAbstractions.@kernel function U1plaquette!(plx, U, Nx, Ny)

    i1, i2 = KernelAbstractions.@index(Global, NTuple)

    iu1 = mod(i1, Nx) + 1
    iu2 = mod(i2, Ny) + 1
    
    plx[i1,i2] = real(U[i1,i2,1] *
                      U[iu1,i2,2] *
                      conj(U[i1,iu2,1] *
                           U[i1,i2,2]))
    
end

# function qtop!(plx, U, prm::LattParm)

#     i1 = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#     i2 = (CUDA.blockIdx().y - 1) * CUDA.blockDim().y + CUDA.threadIdx().y

#     iu1 = mod(i1, prm.iL[1]) + 1
#     iu2 = mod(i2, prm.iL[2]) + 1
    
#     plx[i1,i2] = CUDA.angle(U[i1,i2,1] *
#                             U[iu1,i2,2] *
#                             conj(U[i1,iu2,1] *
#                                  U[i1,i2,2]))
    
#     return nothing
# end

# function force!(U1ws::U1, lp::LattParm)
#     CUDA.@sync begin
#         CUDA.@cuda threads=lp.kprm.threads blocks=lp.kprm.blocks U1force!(U1ws.frc1, U1ws.frc2, U1ws.U, lp.beta, lp.iL[1], lp.iL[2])
#     end
#     return nothing
# end

# function U1force!(frc1, frc2, U, beta, Nx, Ny)
    
#     i1 = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#     i2 = (CUDA.blockIdx().y - 1) * CUDA.blockDim().y + CUDA.threadIdx().y

#     iu1 = mod(i1, Nx) + 1
#     iu2 = mod(i2, Ny) + 1

#     v = beta * imag(U[i1,i2,1] * U[iu1,i2,2] * conj(U[i1,iu2,1] * U[i1,i2,2]))
    
#     frc1[i1,i2,1]  = -v 
#     frc1[i1,i2,2]  =  v 
#     frc2[iu1,i2,2] = -v 
#     frc2[i1,iu2,1] =  v 

#     return nothing
# end

function force!(U1ws::U1, lp::LattParm)
    event = U1force!(lp.device)(U1ws.frc1, U1ws.frc2, U1ws.U, lp.beta, lp.iL[1], lp.iL[2], ndrange=(lp.iL[1], lp.iL[2]), workgroupsize=lp.kprm.threads)
    wait(event)
    return nothing
end

KernelAbstractions.@kernel function U1force!(frc1, frc2, U, beta, Nx, Ny)
    
    i1, i2 = KernelAbstractions.@index(Global, NTuple)

    iu1 = mod(i1, Nx) + 1
    iu2 = mod(i2, Ny) + 1

    v = beta * imag(U[i1,i2,1] * U[iu1,i2,2] * conj(U[i1,iu2,1] * U[i1,i2,2]))
    
    frc1[i1,i2,1]  = -v 
    frc1[i1,i2,2]  =  v 
    frc2[iu1,i2,2] = -v 
    frc2[i1,iu2,1] =  v 
end

# function action(U1ws::U1workspace, lp::LattParm)
#     return U1action(U1ws.U, lp.beta, lp.iL[1], lp.iL[2], lp.kprm.threads, lp.kprm.blocks)
# end

# function U1action(U, beta, Nx, Ny, threads, blocks)
#     plaquettes = CUDA.CuArray{Float64}(undef, Nx, Ny)
#     return U1action(plaquettes, U, beta, Nx, Ny, threads, blocks)
# end

# function U1action(plaquettes, U, beta, Nx, Ny, threads, blocks)
#     CUDA.@sync begin
#         CUDA.@cuda threads=threads blocks=blocks U1plaquette!(plaquettes, U, Nx, Ny)
#     end
#     S = beta * ( Nx * Ny - reduce(+, plaquettes) )

#     return S
# end


function action(U1ws::U1workspace, lp::LattParm)
    return U1action(U1ws.U, lp.beta, lp.iL[1], lp.iL[2], lp.device, lp.kprm.threads, lp.kprm.blocks)
end

function U1action(U, beta, Nx, Ny, device, threads, blocks)
    plaquettes = to_device(device, zeros(Float64, Nx, Ny))
    return U1action(plaquettes, U, beta, Nx, Ny, device, threads, blocks)
end

function U1action(plaquettes, U, beta, Nx, Ny, device, threads, blocks)
    event = U1plaquette!(device)(plaquettes, U, Nx, Ny, ndrange=(Nx, Ny), workgroupsize=threads)
    wait(event)
    S = beta * ( Nx * Ny - reduce(+, plaquettes) )

    return S
end
