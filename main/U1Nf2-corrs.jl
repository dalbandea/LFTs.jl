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

obs = U1PCACCorrelator(lp0)

measure(obs, A0, lp0)
