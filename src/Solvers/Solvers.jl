abstract type AbstractSolver end

struct CG <: AbstractSolver
    maxiter
    tol
    r
    p
    Ap
    tmp
    function CG(maxiter, tol, X)
        r = similar(X)
        p = similar(X)
        Ap = similar(X)
        tmp = similar(X)
        return new(maxiter, tol, r, p, Ap, tmp)
    end
end


invert!(solver::CG, so, A::Function, si, lftws::LFTworkspace, lp::LattParm) = cg!(so, A, si, solver, lftws, lp)


function cg!(so, A::Function, si, solver::CG, lftws::LFTworkspace, lp::LattParm)

    r  = solver.r
    p  = solver.p
    Ap = solver.Ap
    tmp = solver.tmp
    
    so .= zero(eltype(so))
    r  .= si
    p  .= si
    norm = mapreduce(x -> abs2(x), +, si)
    err = zero(lftws.PRC)
    
	    # println( tol)
	    iterations = 0
    for i in 1:solver.maxiter
        A(Ap, tmp, p, lftws, lp)
        prod  = dot(p, Ap)
        alpha = norm/prod

        so .= so .+ alpha .*  p
        r  .= r  .- alpha .* Ap

        err = mapreduce(x -> abs2(x), +, r)
        
        if err < solver.tol
		iterations=i
            break
        end

        beta = err/norm
        p .= r .+ beta .* p
        
        norm = err;
    end

    if err > solver.tol
	    println(err)
        error("CG not converged after $(solver.maxiter) iterationss")
    end
    
    return iterations
end
