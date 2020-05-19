
"""
Tests for initialization portion of algorithm: spherical harmonics, box lists
"""
using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using Random, Distributions
import FMMSteps
import FMM3D
import BoxStruct
import GridStruct

function sphericalHarmonicsTest()
    α = 0.5*1/sqrt(pi)
    β = 0.5*α
    θ = 0.0
    ϕ = 1.0
    s = sin(θ)
    c = cos(θ)
    expo = ℯ^(im*ϕ)
    expt = ℯ^(2*im*ϕ)
    nexpo = ℯ^(-im*ϕ)
    nexpt = ℯ^(-2*im*ϕ)
    Ylma = [α, α*sqrt(3/2)*s*nexpo, α*sqrt(3)*c, -α*sqrt(3/2)*s*expo]
    Ylmb = [β*sqrt(15/2)*s^2*nexpt, β*sqrt(15/2)*s*c*nexpo, β*sqrt(5)*(3c^2-1), -β*sqrt(15/2)*s*c*expo, β*sqrt(15/2)*s^2*expt]
    Ylm_vec_theta_expected_a = convert(Array{ComplexF64,1},Ylma)
    Ylm_vec_theta_a = spherical_harmonics(1,θ,ϕ)

    println("Testing Spherical Harmonics, lmax = 1")
    @test Ylm_vec_theta_a == Ylm_vec_theta_expected_a

    Ylm_vec_theta_expected_b = convert(Array{ComplexF64,1},[Ylma;Ylmb])
    Ylm_vec_theta_b = spherical_harmonics(2,θ,ϕ)
    println("Testing Spherical Harmonics, lmax = 2 ")
    @test Ylm_vec_theta_b == Ylm_vec_theta_expected_b

    Ylma_rev = [α, -α*sqrt(3/2)*s*expo, α*sqrt(3)*c, α*sqrt(3/2)*s*nexpo]
    Ylmb_rev = [β*sqrt(15/2)*s^2*expt, -β*sqrt(15/2)*s*c*expo, β*sqrt(5)*(3c^2-1), β*sqrt(15/2)*s*c*nexpo, β*sqrt(15/2)*s^2*nexpt]
    Ylm_vec_theta_expected_rev_a = convert(Array{ComplexF64,1},Ylma_rev)
    Ylm_vec_theta_rev_a = spherical_harmonics(1,θ,ϕ,true)

    println("Testing Spherical Harmonics reversed ms, lmax = 1")
    @test Ylm_vec_theta_rev_a == Ylm_vec_theta_expected_rev_a

    Ylm_vec_theta_expected_rev_b = convert(Array{ComplexF64,1},[Ylma_rev;Ylmb_rev])
    Ylm_vec_theta_rev_b = spherical_harmonics(2,θ,ϕ,true)
    println("Testing Spherical Harmonics reversed ms, lmax = 2")
    @test Ylm_vec_theta_rev_b == Ylm_vec_theta_expected_rev_b

    println("Testing Spherical Harmonics reversed ms indices")
    ind1 = lm_to_spherical_harmonic_index(1,1)
    negmind1 = lm_to_spherical_harmonic_index(1,-1)
    @test Ylm_vec_theta_expected_a[ind1] == Ylm_vec_theta_expected_rev_a[negmind1]
end

sphericalHarmonicsTest()

function listTests()

    grid = FMM3D.step0(particles,2,4)

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
