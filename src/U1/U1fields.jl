

@doc raw"""
    struct U1workspace{T}

Allocates all the necessary fields for a HMC simulation of a U(1) model:

- `PRC`: precision; must be `ComplexF64`.
- `U`: ``U`` gauge field.
- `frc`
"""
struct U1workspace{T} <: U1
    PRC::Type{T}
    U
    frc1
    frc2
    mom
    function U1workspace(::Type{T}, lp::U1Parm) where {T <: Complex}
        U = CUDA.ones(T, lp.iL..., 2)
        frc1 = CUDA.ones(Float64, lp.iL..., 2)
        frc2 = CUDA.ones(Float64, lp.iL..., 2)
        mom = CUDA.ones(Float64, lp.iL..., 2)
        return new{T}(T, U, frc1, frc2, mom)
    end
end

function copy!(U1ws_dst::U1, U1ws_src::U1, lp::U1Parm)
    U1ws_dst.U .= U1ws_src.U
    return nothing
end

# @doc raw"""
#     function randomize!(phiws, lp)

# Randomizes:

# - `phiws.phi` from Gaussian.
# """
# function randomize!(phiws::Phi4, lp::Phi4Parm)
#     phiws.phi .= Random.randn(phiws.PRC, size(phiws.phi)...)
#     return nothing
# end

