
abstract type Phi4 <: LFTworkspace end

struct Phi4Parm <: LattParm
    iL::Tuple{Int64,Int64}
    beta::Float64
    lambda::Float64
end
export Phi4Parm

include("Phi4fields.jl")
export Phi4workspace, randomize! 

include("Phi4action.jl")
export action, force!

include("Phi4hmc.jl")
include("Phi4checks.jl")
export check_force, reversibility!

include("Phi4measurements.jl")
export magnetization, susceptibility, chi2, correlation_function, correlation_function2, G0, correlation_matrix

include("Phi4gradientflow.jl")
export he_flow

include("Phi4Upscaled/Phi4Upscaled.jl")

# Glossary of variable name meanings

# ws = workspace
# lp = lattice parameter
# frc = force

