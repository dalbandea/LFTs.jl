using LFTs
using LinearAlgebra
using CUDAKernels
import CUDA
CUDA.allowscalar(true)

# U1 Theory parameters
lsize = 12
beta = 5.0
mass = 0.6

model = LFTs.U1Nf2(Float64, iL = (lsize, lsize), beta = beta, am0 = mass)

alg = HMC(
          integrator = Leapfrog(1.0, 10),
          ntherm = 10,
          ntraj = 100,
         )

hmcws = LFTs.sampler(model, alg)

randomize!(model)
LFTs.generate_pseudofermions!(model, hmcws)

@testset "Pseudofermion generation" begin
    xi = similar(hmcws.X)

    invert!(xi, gamm5Dw_sqr_msq!, hmcws.F, model.sws, model)

    @test isapprox(dot(xi, hmcws.F), dot(hmcws.X, hmcws.X))
end


@testset "HMC reversibility" begin
    model_bckp = deepcopy(model)

    reversibility!(model, hmcws)

    ΔU = model.U .- model_bckp.U

    @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, ΔU), atol = 1e-15)
end

get_field(u1ws::LFTs.U1Nf2) = u1ws.U


function analytic_force(u1ws::LFTs.U1Nf2, hmcws::LFTs.U1Nf2HMC)
    force!(u1ws, hmcws)
    frc = hmcws.frc1 .+ hmcws.frc2 .+ hmcws.pfrc
    return frc
end


function infinitesimal_transformation(field_elem, epsilon, lftws2::LFTs.U1)
    return field_elem * exp(im*epsilon)
end


@testset "HMC force" begin
    ΔF = force_test(model, hmcws, 1e-5)

    @test isapprox(zero(model.PRC), ΔF, atol = 1e-5)
end
