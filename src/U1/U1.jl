
abstract type U1 <: LFTworkspace end

import CUDA

struct KernelParm
    threads::Tuple{Int64,Int64}
    blocks::Tuple{Int64,Int64}
end
export KernelParm

struct U1Parm <: LattParm
    iL::Tuple{Int64,Int64}
    beta::Float64
    kprm::KernelParm
end
export U1Parm


include("U1fields.jl")
export U1workspace

include("U1action.jl")
export action, U1plaquette!

include("U1hmc.jl")
# include("Phi4checks.jl")
export Hamiltonian, generate_momenta!, update_fields!, U1_update_field!, update_momenta!


# Glossary of variable name meanings

# ws = workspace
# lp = lattice parameter
# frc = force


