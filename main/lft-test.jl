using Revise
import Pkg
Pkg.activate(".")
using LFTs

# CP Theory parameters
N = 2
lsize = 12
beta = 1.1
lp0 = CPNParm(N, (lsize,lsize), beta)

# Initialize CPN workspace
A0 = CPworkspace(Float64, lp0)

# Initialize fields
randomize!(A0, lp0)

# HMC
tau = 1.0
ns = 2
integrator = OMF4(tau, ns)

@time hmc!(A0, integrator, lp0)

for i in 1:100000
    println(i)
    dH = hmc!(A0, epsilon, ns, lp0)
    S = action(A0, lp0)

    global io_stat = open("test_output.txt", "a")
    write(io_stat, "$(S),$(dH)\n")
    close(io_stat)
end


# Phi4

# Phi4 Theory parameters
lsize = 16
beta = 0.634
lambda = 0.5
lp0 = Phi4Parm((lsize,lsize), beta, lambda)

# Initialize CPN workspace
A0 = Phi4workspace(Float64, lp0)

A0 = Phi4Upscaledworkspace(Float64, lp0)

# Initialize fields
randomize!(A0, lp0)

# HMC
tau = 1.0
ns = 10
integrator = Leapfrog(tau, ns)

@time hmc!(A0, integrator, lp0)

mags = Vector{Float64}()

for i in 1:10000
    println(i)
    dH = hmc!(A0, integrator, lp0)

    # push!(mags, magnetization(A0.phi))
    # S = action(A0, lp0)

    # global io_stat = open("test_output.txt", "a")
    # write(io_stat, "$(S),$(dH)\n")
    # close(io_stat)
end

plot(mags)

sampler = HMC(integrator = integrator)

sample!(A0, sampler, lp0)

sampler.ntherm = 10
sampler.ntraj = 10000
sampler.thermalized = true
sampler.nmeas = 10

sample_and_measure!("nothing", A0, sampler, lp0, verbose = true)


# U1 quenched
using Revise
import Pkg
Pkg.activate(".")
using LFTs
using CUDAKernels

# U1 Theory parameters
lsize = 12
beta = 5.0
kprm = KernelParm((lsize, 1), (1,lsize))
lp0 = U1Parm((lsize, lsize), beta, CUDADevice(), kprm)


A0 = U1quenchedworkspace(ComplexF64, lp0)

@time generate_momenta!(A0, lp0)

@time force!(A0, lp0)

update_fields!(A0, 0.001, lp0)


# HMC
tau = 1.0
ns = 10
integrator = Leapfrog(tau, ns)

@time hmc!(A0, integrator, lp0)


# U1 Nf=2
using Revise
import Pkg
Pkg.activate(".")
using LFTs
using CUDAKernels

# U1 Theory parameters
lsize = 12
beta = 5.0
kprm = KernelParm((lsize, 1), (1,lsize))
lp0 = U1Parm((lsize, lsize), beta, CUDADevice(), kprm)
mass = 0.6

A0 = U1Nf2workspace(Float64, lp0, mass)

@time generate_momenta!(A0, lp0)

@time force!(A0, lp0)

update_fields!(A0, 0.001, lp0)

@time generate_pseudofermions!(A0, lp0)

xi = similar(A0.X)

invert!(A0.sws, xi, gamm5Dw_sqr_msq!, A0.F, A0, lp0)

isapprox(dot(xi, A0.F), dot(A0.X, A0.X))


# HMC
tau = 1.0
ns = 10
integrator = Leapfrog(tau, ns)

@time hmc!(A0, integrator, lp0)
