
function generate_momenta!(U1ws::U1, lp::U1Parm)
    # Create momenta for U1
    U1ws.mom .= to_device(lp.device, randn(Float64, size(U1ws.mom)))
    return nothing
end

function Hamiltonian(U1ws::U1, lp::U1Parm)
    return Hamiltonian(U1ws.mom, U1ws, lp)
end

function Hamiltonian(mom, U1ws::U1, lp::U1Parm)
    H = CUDA.mapreduce(x -> x^2, +, mom)/2.0 + action(U1ws, lp)
    return H
end

function update_momenta!(U1ws::U1, epsilon, lp::U1Parm)
    force!(U1ws, lp)
    update_momenta!(U1ws.mom, U1ws.frc1, U1ws.frc2, epsilon, U1ws::U1)
    return nothing
end

function update_momenta!(mom, frc1, frc2, epsilon, U1ws::U1)
    mom .= mom .+ epsilon * (frc1 .+ frc2)
    return nothing
end


# function update_fields!(U1ws::U1, epsilon, lp::U1Parm)
#     CUDA.@sync begin
#         CUDA.@cuda threads=lp.kprm.threads blocks=lp.kprm.blocks U1_update_field!(U1ws.U, U1ws.mom, epsilon)
#     end
#     return nothing
# end


# function U1_update_field!(U, mom, eps)

#     i1 = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#     i2 = (CUDA.blockIdx().y - 1) * CUDA.blockDim().y + CUDA.threadIdx().y

#     for id in 1:2
#         U[i1,i2,id] = complex(CUDA.cos(eps*mom[i1,i2,id]), CUDA.sin(eps*mom[i1,i2,id])) * U[i1,i2,id]
#     end
    
#     return nothing
# end

function update_fields!(U1ws::U1, epsilon, lp::U1Parm)
    event = U1_update_field!(CUDAKernels.CUDADevice())(U1ws.U, U1ws.mom, epsilon, ndrange=(lp.iL[1], lp.iL[2]), workgroupsize=lp.kprm.threads)
    wait(event)
    return nothing
end


KernelAbstractions.@kernel function U1_update_field!(U, mom, epsilon)

    i1, i2 = KernelAbstractions.@index(Global, NTuple)

    for id in 1:2
        U[i1,i2,id] = complex(CUDA.cos(epsilon*mom[i1,i2,id]), CUDA.sin(epsilon*mom[i1,i2,id])) * U[i1,i2,id]
    end
end
