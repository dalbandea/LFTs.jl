
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




# Glossary of variable name meanings

# ws = workspace
# lp = lattice parameter
# frc = force

