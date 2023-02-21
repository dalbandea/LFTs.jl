using LFTs
using LinearAlgebra

I = 0.5
iT = 64


model = QuantumRotor(I = I, iT = iT)
randomize!(model)


sampler = HMC(
              integrator = Leapfrog(1.0, 10),
              ntherm = 10,
              ntraj = 100,
             )

hmcws = LFTs.sampler(model, sampler)


@testset "HMC reversibility" begin
    model_bckp = deepcopy(model)

    reversibility!(model, hmcws)

    Δϕ = model.phi .- model_bckp.phi

    @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, Δϕ), atol = 1e-15)
end

get_field(qrws::LFTs.QuantumRotor) = qrws.phi


function analytic_force(qrws::LFTs.QuantumRotor, hmcws::LFTs.QuantumRotorHMC)
    force!(qrws, hmcws)
    return hmcws.frc
end


function infinitesimal_transformation(field_elem, epsilon, lftws2::LFTs.QuantumRotor)
    return field_elem + epsilon
end


@testset "HMC force" begin
    ΔF = force_test(model, hmcws, 1e-5)

    @test isapprox(zero(model.PRC), ΔF, atol = 1e-5)
end

