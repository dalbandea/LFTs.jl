abstract type AbstractCorrelator <: AbstractObservable end
abstract type U1AbstractCorrelator <: AbstractCorrelator end

struct U1PionCorrelator <: U1AbstractCorrelator
    name::String
    ID::String
    filepath::String
    R1
    R2
    S
    S0
    C # correlator
    function U1PionCorrelator(lp::U1Parm; wdir::String = "./results/trash/", 
                               name::String = "U(1) pion correlator with Nf=2", 
                               ID::String = "corr", 
                               mesdir::String = "measurements/", 
                               extension::String = ".txt")
        filepath = joinpath(wdir, mesdir, ID*extension)
        R1 = to_device(lp.device, zeros(complex(Float64), lp.iL..., 2))
        R2 = copy(R1)
        S = copy(R1)
        S0 = copy(R1)
        C = zeros(Float64, lp.iL[1])
        return new(name, ID, filepath, R1, R2, S, S0, C)
    end
end
export U1PionCorrelator


function invert_sources!(corrws::U1AbstractCorrelator, U1ws::U1Nf2, lp::U1Parm)

    S0 = corrws.S0
    S = corrws.S
    R1 = corrws.R1
    R2 = corrws.R2

    # Source 1
    S0 .= zero(eltype(S0))
    CUDA.@allowscalar S0[1,1,1] = one(eltype(S0))
    S = similar(S0)

    ## Solve g5D S = S0 for S
	iter = invert!(U1ws.sws, S, gamm5Dw_sqr_msq!, S0, U1ws, lp)
    gamm5Dw!(R1, S, U1ws, lp)

    # Source 2
    S0 .= zero(eltype(S0))
    CUDA.@allowscalar S0[1,1,2] = one(eltype(S0))

	iter = invert!(U1ws.sws, S, gamm5Dw_sqr_msq!, S0, U1ws, lp)
    gamm5Dw!(R2, S, U1ws, lp)

    # NOTE: R1[1,1,1]-R2[1,1,2] is the chiral condensate. Checking that it has
    # no imaginary part is a good test
    
    return nothing
end
export invert_sources!


function pion_correlator_function(corrws::U1AbstractCorrelator, t, lp::U1Parm)

    Ct = zero(ComplexF64)
    a = zero(ComplexF64)
    b = zero(ComplexF64)
    c = zero(ComplexF64)
    d = zero(ComplexF64)

    # NOTE: this should be ultraslow. It may be better to put R1 and R2 into the
    # CPU prior to calling this function. For GPU, the best one can do is to
    # reduce columns of the GPU array.
    for x in 1:lp.iL[1]
        CUDA.@allowscalar begin 
            a = corrws.R1[x,t,1]
            b = corrws.R1[x,t,2]
            c = corrws.R2[x,t,1]
            d = corrws.R2[x,t,2]
        end
        Ct += abs(dot(a,a) + dot(b,b) + dot(c,c) + dot(d,d))
    end

    return Ct
end


function pion_correlator_function(corrws::U1AbstractCorrelator, lp::U1Parm)
    for t in 1:lp.iL[1]
        corrws.C[t] = pion_correlator_function(corrws, t, lp) |> real
    end
end

function measure(corrws::U1PionCorrelator, U1ws::U1Nf2, lp::U1Parm)
    invert_sources!(corrws, U1ws, lp)
    pion_correlator_function(corrws, lp)
    return nothing
end

struct U1PCACCorrelator <: U1AbstractCorrelator
    name::String
    ID::String
    filepath::String
    R1
    R2
    S
    S0
    C # correlator
    function U1PCACCorrelator(lp::U1Parm; wdir::String = "./results/trash/", 
                               name::String = "U(1) PCAC correlator with Nf=2", 
                               ID::String = "corr", 
                               mesdir::String = "measurements/", 
                               extension::String = ".txt")
        filepath = joinpath(wdir, mesdir, ID*extension)
        R1 = to_device(lp.device, zeros(complex(Float64), lp.iL..., 2))
        R2 = copy(R1)
        S = copy(R1)
        S0 = copy(R1)
        C = zeros(Float64, lp.iL[1])
        return new(name, ID, filepath, R1, R2, S, S0, C)
    end
end
export U1PCACCorrelator

function pcac_correlation_function(corrws::U1AbstractCorrelator, t, lp::U1Parm)

    Ct = zero(ComplexF64)
    a = zero(ComplexF64)
    b = zero(ComplexF64)
    c = zero(ComplexF64)
    d = zero(ComplexF64)

    for x in 1:lp.iL[1]
        CUDA.@allowscalar begin 
            a = corrws.R1[x,t,1]
            b = corrws.R1[x,t,2]
            c = corrws.R2[x,t,1]
            d = corrws.R2[x,t,2]
        end
        Ct += -imag(a*conj(c)) - imag(b*conj(d))
    end
    Ct *= 2

    # NOTE: another test would be to check if there is imaginary part of Ct

    return Ct
end

function pcac_correlation_function(corrws::U1AbstractCorrelator, lp::U1Parm)
    for t in 1:lp.iL[1]
        corrws.C[t] = pcac_correlation_function(corrws, t, lp)
    end
end


function measure(corrws::U1AbstractCorrelator, U1ws::U1Nf2, lp::U1Parm)
    invert_sources!(corrws, U1ws, lp)
    correlation_function(corrws, lp)
    return nothing
end

correlation_function(corrws::U1PionCorrelator, lp::U1Parm) = pion_correlator_function(corrws, lp)
correlation_function(corrws::U1PCACCorrelator, lp::U1Parm) = pcac_correlation_function(corrws, lp)
