using Test
using LFTs

include("hmc-test.jl")

@testset verbose = true "U1 Nf=2" begin
    include("U1/U1tests.jl")
end
