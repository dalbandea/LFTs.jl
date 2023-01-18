using Revise, Dates
import Pkg
Pkg.activate(".")
using LFTs

lsize = 16
beta = 0.634
lambda = 0.5
lp0 = Phi4Parm((lsize,lsize), beta, lambda)

# Initialize CPN workspace
A0 = Phi4workspace(Float64, lp0)

# Initialize fields
randomize!(A0, lp0)

# HMC
tau = 1.0
ns = 15
integrator = Leapfrog(tau, ns)
ntherm = 1000
ntraj = 10000
nmeas = 20
sampler = HMC(integrator=integrator, ntherm=ntherm, ntraj=ntraj, nmeas=nmeas)

# Create dirs
wdir = joinpath("./results/trash/", string(now()))
mesdir = joinpath(wdir, "measurements")
anadir = joinpath(wdir, "analysis")
mkpath(wdir)
mkpath(mesdir)

# Logs for reproducibility
reproducibility_log(wdir)

# Define observables
observables = [
               Phi4AverageOnPoint(wdir=wdir),
              ]

# Sample and measure, outputing to file
sample_and_measure!(observables, A0, sampler, lp0)

# Analyze
analyze.(observables, wdir=anadir)
