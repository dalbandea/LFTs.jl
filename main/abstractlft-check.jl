using Revise
import Pkg
Pkg.activate(".")
using LFTs
using CUDAKernels

# Phi4

beta = 0.576
lambda = 0.5
lsize = 8

model = Phi4(
             beta=beta, 
             lambda=lambda, 
             iL=(8,8)
            )

sampler = HMC(
              integrator = Leapfrog(
                                   1.0,
                                   10
                                  ),
              ntherm = 10,
              ntraj = 100,
             )

@code_warntype sample!(model, samplerws)


# U1

beta = 0.2
lsize = 8

model = LFTs.U1Quenched(Float64, beta = beta, iL = (lsize, lsize))


sampler = HMC(
              integrator = Leapfrog(
                                   1.0,
                                   10
                                  ),
              ntherm = 10,
              ntraj = 100,
             )

samplerws = LFTs.sampler(model, sampler)

sample!(model, samplerws)


# U1Nf2

beta = 0.2
lsize = 8
mass = 0.6

model = LFTs.U1Nf2(Float64, beta = beta, iL = (lsize, lsize), am0 = mass)

sampler = HMC(
              integrator = Leapfrog(1.0, 10),
              ntherm = 10,
              ntraj = 100,
             )

samplerws = LFTs.sampler(model, sampler)

for i in 1:10
    sample!(model, samplerws)
end

bck = copy(model.U)

LFTs.molecular_dynamics!(model, samplerws)
samplerws.mom .= .-samplerws.mom

LFTs.molecular_dynamics!(model, samplerws)

bck .- model.U



