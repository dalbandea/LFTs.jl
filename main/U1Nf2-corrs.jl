# U1 Nf=2
using Revise
import Pkg
Pkg.activate(".")
using LFTs
using CUDAKernels
using Dates

# U1 Theory parameters
lsize = 16
beta = 2.0
kprm = KernelParm((lsize, 1), (1,lsize))
lp0 = U1Parm((lsize, lsize), beta, CUDADevice(), kprm)
mass = -0.18840579710144945

A0 = U1Nf2workspace(Float64, lp0, mass)

# HMC
tau = 1.0
ns = 30
integrator = Leapfrog(tau, ns)
ntherm = 100
ntraj = 10000
nmeas = 1
sampler = HMC(integrator=integrator, ntherm=ntherm, ntraj=ntraj, nmeas=nmeas)

# Create dirs
# wdir = joinpath("/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/1-Phi4/2-distribution-between-neighbors", string(now()))
wdir = joinpath("./results/trash/U1CriticalMassCheck", string(now()))
mesdir = joinpath(wdir, "measurements")
anadir = joinpath(wdir, "analysis")
mkpath(wdir)
mkpath(mesdir)
mkpath(anadir)

# Logs for reproducibility
reproducibility_log(wdir)

# Define observables
observables = [
               U1PionCorrelator(lp0, wdir=wdir),
               U1PCACCorrelator(lp0, wdir=wdir),
              ]

# Sample and measure, outputing to file
sample_and_measure!(observables, A0, sampler, lp0, verbose = true, get_history = false)

# Analyze
# analyze.(observables, wdir=anadir)
