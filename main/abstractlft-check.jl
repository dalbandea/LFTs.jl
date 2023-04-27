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
samplerws = LFTs.sampler(model, sampler)

randomize!(model)

@code_warntype sample!(model, samplerws)

dHs = Vector{Float64}()

for i in 1:300000
    dH = hmc!(model, samplerws)
    push!(dHs, dH)
end

using ADerrors

phitest = "testphi14"
uwdh = uwreal(exp.(-dHs)[10000:3:end], phitest)
uwerr(uwdh)
uwdh

uwM = uwreal(abs.(Ms)[10000:end], phitest)
uwerr(uwM)
uwM

taui(uwdh ,phitest)
dtaui(uwdh ,phitest)


# U1

beta = 5.0
lsize = 12
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

# Qs = []
Ss = []


for i in 1:1000000
    println(i)
    @time sample!(model, samplerws)
    # push!(Qs, LFTs.top_charge(model))
    # push!(Ss, LFTs.action(model,samplerws))
    global io = open("HMC-Qs.txt", "a")
    println(io, LFTs.top_charge(model))
    close(io)
end


uwchi = uwreal(Qs.^2/12^2, "test")
uwerr(uwchi)
uwchi
taui(uwchi, "test")
dtaui(uwchi, "test")

using ADerrors

uws = uwreal(Ss.^1, "test2")
uwerr(uws)
uws
taui(uws, "test2")
dtaui(uws, "test2")

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
