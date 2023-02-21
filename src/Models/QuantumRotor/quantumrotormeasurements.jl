
abstract type Susceptibility <: AbstractObservable end

struct QRMasterObs end

mutable struct QRTopologicalCharge <: AbstractScalar
    name::String
    ID::String
    filepath::String
    result::Float64
    history::Vector{Float64}
    function QRTopologicalCharge(; wdir::String = "./results/trash/", 
                              name::String = "Topological charge", 
                              ID::String = "topcharge", 
                              mesdir::String = "measurements/", 
                              extension::String = ".txt")
        filepath = joinpath(wdir, mesdir, ID*extension)
        result = zero(Float64)
        history = Vector{Float64}()
        return new(name, ID, filepath, result, history)
    end
end
export QRTopologicalCharge

function top_charge(qrws::QuantumRotor)
    Q = 0.0

    for i in 1:qrws.params.iT-1
        Q += Mod(qrws.phi[i+1]-qrws.phi[i], 2pi)
    end

    Q += Mod(qrws.phi[1]-qrws.phi[end], 2pi)

    return Q/2pi
end

function (obs::QRTopologicalCharge)(qrws::QuantumRotor)
    obs.result = top_charge(qrws)
    return nothing
end
