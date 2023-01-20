
function magnetization(phi)
    return sum(phi) / length(phi)
end

function G0(phi)
    return sum(phi.^2) / length(phi)
end

function susceptibility(phi)
    return sum(phi.^2) / length(phi)
end

function chi2(mags, lp)
    V = lp.iL[1] * lp.iL[2]
    M = sum(mags) / length(mags)
    M2 = mags.^2
    chi = V * (M2 .- M^2)
    return chi
end

function correlation_function(phi, shift::Tuple{Int64, Int64})
    L = size(phi)[1]
    C = phi .* circshift(phi, shift)

    return sum(C)/L^2
end

circindex(i::Int,N::Int) = 1 + mod(i-1,N) # for symmetric BC

# function correlation_function(phi, y)
#     C = 0.0
#     L = size(phi)[1]

#     for i in 1:L
#         for j in 1:L
#             C += phi[i,j] * (phi[circindex(i+y,L), j] + phi[i, circindex(j+y, L)])
#         end
#     end

#     return C/L^2/2
# end

function correlation_function(phi)
    L = size(phi)[1]
    C = zeros(L)

    for i in 1:L
        for j in 1:L
            C[i] += correlation_function(phi, (i-1, j-1))
        end
    end

    return C ./ L
end

function correlation_matrix(phi)
    L = size(phi)[1]
    C = zeros((L,L))

    for i in 1:L
        for j in 1:L
            C[i,j] = correlation_function(phi, (i-1, j-1))
        end
    end

    return C
end


mutable struct Phi4AverageOnPoint <: AbstractScalar
    name::String
    ID::String
    filepath::String
    result::Float64
    history::Vector{Float64}
    function Phi4AverageOnPoint(; wdir::String = "./results/trash/", 
                              name::String = "Average on point", 
                              ID::String = "avgpt", 
                              mesdir::String = "measurements/", 
                              extension::String = ".txt")
        filepath = joinpath(wdir, mesdir, ID*extension)
        result = zero(Float64)
        history = Vector{Float64}()
        return new(name, ID, filepath, result, history)
    end
end
export Phi4AverageOnPoint

function average_on_point(phiws::Phi4, i::Int64, j::Int64, lp::Phi4Parm)
    iu = mod1(i+1, lp.iL[1])
    id = mod1(i-1, lp.iL[1])
    ju = mod1(i+1, lp.iL[2])
    jd = mod1(i-1, lp.iL[2])

    avg = (phiws.phi[iu,ju] + phiws.phi[id,jd] + phiws.phi[id,ju] + phiws.phi[iu,jd])/4
    # avg = (phiws.phi[iu,j] + phiws.phi[id,j])/2
    val = phiws.phi[i,j]

    return avg-val
end

function (obs::Phi4AverageOnPoint)(phiws::Phi4, lp::Phi4Parm)
    return average_on_point(phiws, 1, 1, lp)
end

function analyze(obs::Phi4AverageOnPoint; wdir::String = "./results/trash/")
    filename = obs.ID*".png"
    data = read(obs)
    pl = Plots.histogram(data)
    Plots.savefig(pl, joinpath(wdir, filename))
    return nothing
end



