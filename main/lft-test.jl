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
ns = 15
integrator = OMF4(tau, ns)

@time HMC!(A0, integrator, lp0)

for i in 1:100000
    println(i)
    dH = HMC!(A0, epsilon, ns, lp0)
    S = action(A0, lp0)

    global io_stat = open("test_output.txt", "a")
    write(io_stat, "$(S),$(dH)\n")
    close(io_stat)
end


# Phi4

# Phi4 Theory parameters
lsize = 8
beta = 0.800
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

@time HMC!(A0, integrator, lp0)

mags = Vector{Float64}()

for i in 1:100000
    println(i)
    dH = HMC!(A0, integrator, lp0)

    push!(mags, magnetization(A0.phi))
    # S = action(A0, lp0)

    # global io_stat = open("test_output.txt", "a")
    # write(io_stat, "$(S),$(dH)\n")
    # close(io_stat)
end

plot(mags)


# U1
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


A0 = U1workspace(ComplexF64, lp0)

generate_momenta!(A0, lp0)

force!(A0, lp0)

update_fields!(A0, 0.001, lp0)


# HMC
tau = 1.0
ns = 10
integrator = Leapfrog(tau, ns)

@time HMC!(A0, integrator, lp0)
