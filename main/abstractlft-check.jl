using Revise
import Pkg
Pkg.activate(".")
using LFTs
using CUDAKernels
using KernelAbstractions

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

beta = 11.25
lsize = 8
device = KernelAbstractions.CPU()

model = LFTs.U1Quenched(Float64, beta = beta, iL = (lsize, lsize), device = device)

# Ucp = copy(model.U)
# model.U .= Ucp

sampler = HMC(
              integrator = Leapfrog(1.0, 10),
              ntherm = 10,
              ntraj = 100,
             )
samplerws = LFTs.sampler(model, sampler)

@time sample!(model, samplerws)

Qs = []
Ss = []

for i in 1:1000
    @time sample!(model, samplerws)
    push!(Qs, LFTs.top_charge(model))
    push!(Ss, LFTs.action(model,samplerws))
end

model_bckp = deepcopy(model)

reversibility!(model, samplerws)

ΔU = model.U .- model_bckp.U

isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, ΔU), atol = 1e-15)

mapreduce(x -> abs2(x), +, ΔU)

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
