using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using EnsembleTests
import FMMSteps
import FMM3D
import BoxStruct
import GridStruct

function sphericalHarmonicsTest()
    θ = pi/2
    ϕ = 1.0
    s = sin(θ)
    c = cos(θ)
    expo = ℯ^(im*ϕ)
    expt = ℯ^(2*im*ϕ)
    nexpo = ℯ^(-im*ϕ)
    nexpt = ℯ^(-2*im*ϕ)
    Ylma = [1.0, sqrt(1/2)*s*nexpo, c, sqrt(1/2)*s*expo]
    Ylmb = [(sqrt(6)/4)*s^2*nexpt , -sqrt(6)/2*s*c*nexpo , (1/2)*(3c^2-1), sqrt(6)/2*s*c*expo, -sqrt(6)/4*s^2*nexpt]
    Ylm_vec_expected_a = convert(Array{ComplexF64,1},Ylma)
    Ylm_vec_a = spherical_harmonics(1,θ,ϕ)

    println("Testing Spherical Harmonics, lmax = 1")
    @test Ylm_vec_a == Ylm_vec_expected_a

    Ylm_vec_expected_b = convert(Array{ComplexF64,1},[Ylma;Ylmb])
    Ylm_vec_b = spherical_harmonics(2,θ,ϕ)
    println("Testing Spherical Harmonics, lmax = 2 ")
    @test Ylm_vec_b == Ylm_vec_expected_b
end

sphericalHarmonicsTest()

function listTests()

    particles, _, minbound, maxbound = particledist1()

    grid = FMM3D.step0(particles,minbound, maxbound,2,4)

    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
    end

    println("Testing L1")
    @test 2 in grid.boxes[2].L1 && 3 in grid.boxes[2].L1 &&  length(grid.boxes[2].L1) == 2
    @test 3 in grid.boxes[3].L1 && 2 in grid.boxes[3].L1 &&  length(grid.boxes[3].L1) == 2
    @test 4 in grid.boxes[4].L1 && 8 in grid.boxes[4].L1 && 6 in grid.boxes[4].L1 &&
        length(grid.boxes[4].L1) == 3
    @test 7 in grid.boxes[7].L1 && 6 in grid.boxes[7].L1 && 8 in grid.boxes[7].L1 &&
        length(grid.boxes[7].L1) == 3
    @test 6 in grid.boxes[6].L1 && 7 in grid.boxes[6].L1 &&  4 in grid.boxes[6].L1 &&
        8 in grid.boxes[6].L1 && length(grid.boxes[6].L1) == 4
    @test 8 in grid.boxes[8].L1 && 7 in grid.boxes[8].L1 &&  6 in grid.boxes[8].L1 &&
        4 in grid.boxes[8].L1 && length(grid.boxes[6].L1) == 4
    @test length(grid.boxes[0].L1) == 0
    @test length(grid.boxes[1].L1) == 0
    @test length(grid.boxes[5].L1) == 0

    println("Testing L2")
    @test 6 in grid.boxes[2].L2 && 7 in grid.boxes[2].L2 &&  length(grid.boxes[2].L2) == 2
    @test 6 in grid.boxes[3].L2 && 7 in grid.boxes[3].L2 &&  length(grid.boxes[3].L2) == 2
    @test 2 in grid.boxes[7].L2 && 3 in grid.boxes[7].L2 &&  length(grid.boxes[7].L2) == 2
    @test 2 in grid.boxes[6].L2 && 3 in grid.boxes[6].L2 &&  length(grid.boxes[6].L2) == 2
    @test length(grid.boxes[0].L2) == 0
    @test length(grid.boxes[1].L2) == 0
    @test length(grid.boxes[4].L2) == 0
    @test length(grid.boxes[8].L2) == 0
    @test length(grid.boxes[5].L2) == 0

    println("Testing L3")
    @test 2 in grid.boxes[8].L3 && 3 in grid.boxes[8].L3 &&  length(grid.boxes[8].L3) == 2
    @test 2 in grid.boxes[4].L3 && 3 in grid.boxes[4].L3 &&  7 in grid.boxes[4].L3 &&
        length(grid.boxes[4].L3) == 3
    @test length(grid.boxes[1].L3) == 0
    @test length(grid.boxes[2].L3) == 0
    @test length(grid.boxes[3].L3) == 0
    @test length(grid.boxes[5].L3) == 0
    @test length(grid.boxes[6].L3) == 0
    @test length(grid.boxes[7].L3) == 0

    println("Testing L4")
    @test 8 in grid.boxes[2].L4 && 4 in grid.boxes[2].L4 && length(grid.boxes[2].L4) == 2
    @test 8 in grid.boxes[3].L4 && 4 in grid.boxes[3].L4 && length(grid.boxes[3].L4) == 2
    @test 4 in grid.boxes[7].L4 && length(grid.boxes[7].L4) == 1
    @test length(grid.boxes[1].L4) == 0
    @test length(grid.boxes[4].L4) == 0
    @test length(grid.boxes[5].L4) == 0
    @test length(grid.boxes[6].L4) == 0
    @test length(grid.boxes[8].L4) == 0
end

listTests()
