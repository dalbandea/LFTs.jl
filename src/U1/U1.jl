
import KernelAbstractions
import CUDA, CUDAKernels
import AMDGPU, ROCKernels

struct KernelParm
    threads::Tuple{Int64,Int64}
    blocks::Tuple{Int64,Int64}
end
export KernelParm

struct U1Parm <: LattParm
    iL::Tuple{Int64,Int64}
    beta::Float64
    device::Union{KernelAbstractions.Device, ROCKernels.ROCDevice}
    kprm::KernelParm
end
export U1Parm


include("U1fields.jl")
export U1quenchedworkspace, U1Nf2workspace

include("U1action.jl")
export action, U1plaquette!, U1action, gauge_action

include("U1hmc.jl")
# include("Phi4checks.jl")
export Hamiltonian, generate_momenta!, update_fields!, U1_update_field!, update_momenta!, generate_pseudofermions!

include("U1dirac.jl")
export U1gamm5Dw!, gamm5Dw!, gamm5Dw_sqr_msq!

to_device(::CUDAKernels.CUDADevice, x) = CUDA.CuArray(x)
to_device(::ROCKernels.ROCDevice, x) = AMDGPU.ROCArray(x)


# Glossary of variable name meanings

# ws = workspace
# lp = lattice parameter
# frc = force


