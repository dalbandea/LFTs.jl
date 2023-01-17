
struct CPNParm <: LattParm
    N::Int64
    iL::Tuple{Int64,Int64}
    beta::Float64
end
export CPNParm

include("CPNfields.jl")
export CPworkspace, randomize!, fill_Lambda!, fill_Jn!, sync_fields!, project_to_Sn!

include("CPNaction.jl")
export action, gauge_frc!, x_frc!, load_frcs!

include("CPNhmc.jl")




# Glossary of variable name meanings

# ws = workspace
# lp = lattice parameter
# frc = force

