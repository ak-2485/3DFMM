
using Test
using CoefHelpers

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
println("Testing Spherical Harmonics")
@test Ylm_vec_theta_a == Ylm_vec_theta_expected_a
Ylm_vec_theta_expected_b = convert(Array{ComplexF64,1},[Ylma;Ylmb])
Ylm_vec_theta_b = spherical_harmonics(2,θ,ϕ)
println("Testing Spherical Harmonics")
@test Ylm_vec_theta_b == Ylm_vec_theta_expected_b

Ylma_rev = [α, -α*sqrt(3/2)*s*expo, α*sqrt(3)*c, α*sqrt(3/2)*s*nexpo]
Ylmb_rev = [β*sqrt(15/2)*s^2*expt, -β*sqrt(15/2)*s*c*expo, β*sqrt(5)*(3c^2-1), β*sqrt(15/2)*s*c*nexpo, β*sqrt(15/2)*s^2*nexpt]
Ylm_vec_theta_expected_rev_a = convert(Array{ComplexF64,1},Ylma_rev)
Ylm_vec_theta_rev_a = spherical_harmonics(1,θ,ϕ,true)
println("Testing Spherical Harmonics reversed")
@test Ylm_vec_theta_rev_a == Ylm_vec_theta_expected_rev_a
Ylm_vec_theta_expected_rev_b = convert(Array{ComplexF64,1},[Ylma_rev;Ylmb_rev])
Ylm_vec_theta_rev_b = spherical_harmonics(2,θ,ϕ,true)
println("Testing Spherical Harmonics reversed")
@test Ylm_vec_theta_rev_b == Ylm_vec_theta_expected_rev_b
println("Testing Spherical Harmonics reversed")
ind1 = lm_to_spherical_harmonic_index(1,1)
negmind1 = lm_to_spherical_harmonic_index(1,-1)
@test Ylm_vec_theta_expected_a[ind1] == Ylm_vec_theta_expected_rev_a[negmind1]
