
# Returns so = γ₅D si in-place.
KernelAbstractions.@kernel function U1gamm5Dw!(so, U, si, am0, Nx, Ny)

    i1, i2 = @index(Global, NTuple)

    iu1 = mod(i1, Nx) + 1
    iu2 = mod(i2, Ny) + 1

    id1 = mod1(i1-1, Nx)
    id2 = mod1(i2-1, Ny)

    A = 0.5 * ( U[i1,i2,1] * (si[iu1,i2,1] - si[iu1,i2,2]) +
                conj(U[id1,i2,1])*(si[id1,i2,1] + si[id1,i2,2]) )
    B = 0.5 * ( U[i1,i2,2] * (si[i1,iu2,1] +
                              complex(-imag(si[i1,iu2,2]),real(si[i1,iu2,2]))) +
                conj(U[i1,id2,2])*(si[i1,id2,1] +
                                    complex(imag(si[i1,id2,2]),-real(si[i1,id2,2]))) )
    A2 = 0.5 * ( U[i1,i2,1] * (si[iu1,i2,1] - si[iu1,i2,2]) -
                conj(U[id1,i2,1])*(si[id1,i2,1] + si[id1,i2,2]) )
    B2 = 0.5 * ( U[i1,i2,2] * (-si[i1,iu2,2] +
                              complex(-imag(si[i1,iu2,1]),real(si[i1,iu2,1]))) -
                conj(U[i1,id2,2])*(si[i1,id2,2] +
                                    complex(-imag(si[i1,id2,1]),real(si[i1,id2,1]))) )

    
    so[i1,i2,1] =  (2.0+am0) * si[i1,i2,1] - A - B
    so[i1,i2,2] = -(2.0+am0) * si[i1,i2,2] - A2 - B2
    
end

function gamm5Dw!(so, si, U1ws::U1Nf2, lp::U1Parm)
    event = U1gamm5Dw!(lp.device)(so, U1ws.U, si, U1ws.am0, lp.iL[1],
                                  lp.iL[2], ndrange=(lp.iL[1], lp.iL[2]),
                                  workgroupsize=lp.kprm.threads)
    wait(event)
    return nothing
end


function gamm5Dw_sqr_msq!(so, tmp, si, U1ws::U1Nf2, lp::U1Parm)

    gamm5Dw!(so, si, U1ws, lp)
    tmp .= so
    gamm5Dw!(so, tmp, U1ws, lp)
    # so .= so .+ (am0^2)
    
    return nothing
end


# function gamm5Dw_sqr(so, U, si, am0::Float64, prm::LattParm, kprm::KernelParm)
#     tmp = similar(so)
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, si, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
    
#     return nothing
# end

# function gamm5Dw_sqr_sqr(so, U, si, am0::Float64, prm::LattParm, kprm::KernelParm)
#     tmp = similar(so)
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, si, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
    
#     return nothing
# end

# function gamm5Dw_sqr_musq(so, tmp, U, si, am0::Float64, mu_j::Float64, prm::LattParm, kprm::KernelParm)

#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, si, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     so .= so .+ (mu_j)^2*si
    
#     return nothing
# end

# function gamm5Dw_sqr_sqr_musq(so, tmp, U, si, am0::Float64, mu_j::Float64, prm::LattParm, kprm::KernelParm)

#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, si, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, tmp, am0, prm)
#     end
#     so .= so .+ (mu_j)^2*si
    
#     return nothing
# end


# """
#    (kernel) gamm5(so, si, am0, prm::LattParm)

# Apply γ₅ to a fermion field. `si` is the input fermion field, and `so` the output.

# # Examples
# ```jldoctest
# julia> CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5(si, so, am0, prm)
# ```
# """
# function gamm5(so, si, am0, prm::LattParm)

#     i1 = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
#     i2 = (CUDA.blockIdx().y - 1) * CUDA.blockDim().y + CUDA.threadIdx().y

#     iu1 = mod(i1, prm.iL[1]) + 1
#     iu2 = mod(i2, prm.iL[2]) + 1

#     id1 = mod1(i1-1, prm.iL[1])
#     id2 = mod1(i2-1, prm.iL[2])

#     so[i1,i2,1] =  si[i1,i2,1]
#     so[i1,i2,2] =  -si[i1,i2,2]
    

#     return nothing
# end


# """
# Dw(so, U, si, am0, prm, kprm)

# Applies γ₅(γ₅D) = D to a fermion field `si` and writes it to `so`.

# # Examples
# julia> Dw(so, U, si, am0, prm, kprm)
# """
# function Dw(so, U, si, am0::Float64, prm::LattParm, kprm::KernelParm)

#     tmp = similar(si)
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5Dw(so, U, si, am0, prm)
#     end
#     tmp .= so
#     CUDA.@sync begin
#         CUDA.@cuda threads=kprm.threads blocks=kprm.blocks gamm5(so, tmp, am0, prm)
#     end
    
#     return nothing
# end


function pf_force!(U1ws::U1Nf2, lp::U1Parm)
    event = U1_tr_dQwdU!(lp.device)(U1ws.pfrc, U1ws.U, U1ws.X, U1ws.g5DX,
                                    lp.iL[1], lp.iL[2], 
                                    ndrange=(lp.iL[1], lp.iL[2]),
                                    workgroupsize=lp.kprm.threads)
    wait(event)
    return nothing
end


KernelAbstractions.@kernel function U1_tr_dQwdU!(frc, U, X, g5DwX, Nx, Ny)

    i1, i2 = @index(Global, NTuple)

    iu1 = mod(i1, Nx) + 1
    iu2 = mod(i2, Ny) + 1

    z1 = conj(U[i1,i2,1])*( conj(X[iu1,i2,1])*(g5DwX[i1,i2,1] + g5DwX[i1,i2,2])  -
                           conj(X[iu1,i2,2])*(g5DwX[i1,i2,1] + g5DwX[i1,i2,2]) ) -
     U[i1,i2,1]*( conj(X[i1,i2,1])*(g5DwX[iu1,i2,1] - g5DwX[iu1,i2,2]) +
                   conj(X[i1,i2,2])*(g5DwX[iu1,i2,1] - g5DwX[iu1,i2,2]) )




    z2 = conj(U[i1,i2,2])*( conj(X[i1,iu2,1])*(g5DwX[i1,i2,1] - (1im)*g5DwX[i1,i2,2])  -
                           conj(X[i1,iu2,2])*((1im)*g5DwX[i1,i2,1] + g5DwX[i1,i2,2])) -
     U[i1,i2,2]*( conj(X[i1,i2,1])*(g5DwX[i1,iu2,1] + (1im)*g5DwX[i1,iu2,2]) +
                   conj(X[i1,i2,2])*((1im)*g5DwX[i1,iu2,1] - g5DwX[i1,iu2,2]) )
   
                   
    
    frc[i1,i2,1] = real((1im)*z1)
    frc[i1,i2,2] = real((1im)*z2)

end
