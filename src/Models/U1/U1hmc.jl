
function generate_momenta!(U1ws::U1, hmcws::AbstractHMC)
    # Create momenta for U1
    hmcws.mom .= to_device(U1ws.params.device, randn(U1ws.PRC, size(hmcws.mom)))
    return nothing
end

function Hamiltonian(U1ws::U1, hmcws::AbstractHMC)
    H = CUDA.mapreduce(x -> x^2, +, hmcws.mom)/2.0 + action(U1ws)
    return H
end

function update_momenta!(U1ws::U1Quenched, epsilon, hmcws::AbstractHMC)
    force!(U1ws, hmcws)
    hmcws.mom .= hmcws.mom .+ epsilon * (hmcws.frc1 .+ hmcws.frc2)
    return nothing
end

# function generate_pseudofermions!(U1ws::U1Nf2, lp::U1Parm)
#     U1ws.X .= to_device(lp.device, randn(complex(U1ws.PRC), lp.iL[1], lp.iL[2], 2))
#     gamm5Dw!(U1ws.F, U1ws.X, U1ws, lp)
#     U1ws.g5DX .= to_device(lp.device, zeros(complex(U1ws.PRC), lp.iL[1], lp.iL[2], 2))
#     return nothing
# end

# function update_momenta!(U1ws::U1Nf2, epsilon, lp::U1Parm)

# 	# Solve DX = F for X
#     iter = invert!(U1ws.sws, U1ws.X, gamm5Dw_sqr_msq!, U1ws.F, U1ws, lp)

# 	# Apply gamm5D to X
#     gamm5Dw!(U1ws.g5DX, U1ws.X, U1ws, lp)
	
# 	# Get fermion part of the force in U1ws.pfrc
#     pf_force!(U1ws, lp)

# 	# Get gauge part of the force in U1ws.frc1 and U1ws.frc2
#     gauge_force!(U1ws, lp)

# 	# Final force is frc1+frc2+frc
#     U1ws.mom .= U1ws.mom .+ epsilon .* (U1ws.frc1 .+ U1ws.frc2 .+ U1ws.pfrc)

# 	return nothing
# end


function update_fields!(U1ws::U1, hmcws::AbstractHMC)
    lp = U1ws.params.lp
    event = U1_update_field!(lp.device)(U1ws.U, hmcws.mom, epsilon, ndrange=(lp.iL[1], lp.iL[2]), workgroupsize=lp.kprm.threads)
    wait(event)
    return nothing
end


KernelAbstractions.@kernel function U1_update_field!(U, mom, epsilon)

    i1, i2 = @index(Global, NTuple)

    for id in 1:2
        U[i1,i2,id] = complex(cos(epsilon*mom[i1,i2,id]), sin(epsilon*mom[i1,i2,id])) * U[i1,i2,id]
    end
end
