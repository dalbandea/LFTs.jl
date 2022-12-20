
@doc raw"""
    struct Phi4Upscaledworkspace{T}

Allocates all the necessary fields for a HMC simulation of a Phi4 model:

- `PRC`: precision; must be `Float16`, `Float32` or `Float64`.
- `phi`: ``\phi^4`` field.
- `frc`
"""
struct Phi4Upscaledworkspace{T, N} <: Phi4
    PRC::Type{T}
    phi::Array{T, N}
    frc::Array{T, N}
    mom::Array{T, N}
    function Phi4Upscaledworkspace(::Type{T}, lp::Phi4Parm) where {T <: AbstractFloat}
        phi = Array{T, 2}(undef, lp.iL...)
        frc = similar(phi)
        mom = similar(phi)
        return new{T, 2}(T, phi, frc, mom)
    end
end

function interpolation_upscale(phi)
    Nx = size(phi, 1)
    Ny = size(phi, 2)
    phi_interp = zeros(eltype(phi), 2*Nx, 2*Ny)

    phi_interp[1:2:end, 1:2:end] .= phi
    phi_interp[1:2:end, 2:2:end] .= 1/2 * (phi .+ circshift(phi,(0,-1)))
    phi_interp[2:2:end, 1:2:end] .= 1/2 * (phi .+ circshift(phi,(-1,0)))
    phi_interp[2:2:end, 2:2:end] .= 1/4 * (phi .+ circshift(phi,(-1,0)) + circshift(phi,(0,-1)) + circshift(phi,(-1,-1)))

    return phi_interp
end

    # upsamples = torch.zeros(sample_shape + [2*self.lsize,2*self.lsize])
    # upsamples[..., ::2, ::2] = samples[...,:,:]
    # upsamples[..., ::2, 1::2] = 1.0/2.0 * (samples + samples.roll(-1,-1))
    # upsamples[..., 1::2, ::2] = 1.0/2.0 * (samples + samples.roll(-1,-2))
    # upsamples[..., 1::2, 1::2] = 1.0/4.0 * (samples + samples.roll(-1,-2) +
    #         samples.roll(-1,-1) + samples.roll((-1,-1),(-2,-1)) )
