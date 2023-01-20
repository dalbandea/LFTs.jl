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

# HMC
tau = 1.0
ns = 10
integrator = Leapfrog(tau, ns)

@time hmc!(A0, integrator, lp0)

# obs = U1PCACCorrelator(lp0)
obs = U1PionCorrelator(lp0)

measure(obs, A0, lp0)


# HMC
tau = 1.0
ns = 15
integrator = Leapfrog(tau, ns)
ntherm = 10
ntraj = 10
nmeas = 1
sampler = HMC(integrator=integrator, ntherm=ntherm, ntraj=ntraj, nmeas=nmeas)

# Create dirs
# wdir = joinpath("/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/1-Phi4/2-distribution-between-neighbors", string(now()))
wdir = joinpath("./results/trash/", string(now()))
mesdir = joinpath(wdir, "measurements")
anadir = joinpath(wdir, "analysis")
mkpath(wdir)
mkpath(mesdir)
mkpath(anadir)

# Logs for reproducibility
reproducibility_log(wdir)

# Define observables
observables = [
               U1PCACCorrelator(lp0, wdir=wdir),
              ]

# Sample and measure, outputing to file
sample_and_measure!(observables, A0, sampler, lp0, verbose = true, get_history = false)

# Analyze
analyze.(observables, wdir=anadir)
