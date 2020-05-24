"""
Helper functions for bulding the 3D FMM coefficients. Including spherical harmonics
which are not my own
(see https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/src/specialfunctions.jl)
and the inner and outer functions for translation and conversion operators.
"""
module CoefHelpers

using GSL

export spherical_harmonics, spherical_harmonic_indices
export lm_to_spherical_harmonic_index, anms, Inm, Onm

function spherical_harmonic_indices(lmax::Int64)
    """
    Returns the arrays for indexing spherical harmonics. This work is not my own.
    (see https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/src/specialfunctions.jl)
    """
    ls = [l for l in 0:lmax for m in -l:l]
    ms = [m for l in 0:lmax for m in -l:l]

    return ls, ms
end

function associated_legendre_indices(lmax::Int64)
    """
    Returns the arrays for indexing the associated legendre polynomials.
    This work is not my own.
    (see https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/src/specialfunctions.jl)
    """
    ls = [l for l in 0:lmax for m in 0:l]
    ms = [m for l in 0:lmax for m in 0:l]

    return ls, ms
end

lm_to_spherical_harmonic_index(l::Int64,m::Int64)::Int64 = l^2 + m + l + 1

function spherical_harmonics(l_max::Int64, θ::Float64, φ::Float64)
    """
    Returns a vector of all spherical harmonics with degree `l <= l_max` as
    used in the 3D FMM; pay attention to normalization for the FMM being
    implemented as these differ in the literature.
    The degree and order (indices) of the elements of the vector are given by
    `spherical_harmonics_indices(l_max::Int)`.
    The associated legendre polynomials are taken from the package GSL.jl.
    This work is derived from
    (see https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/src/specialfunctions.jl)
    and in essence is not my own.
    """

    ls, ms = associated_legendre_indices(l_max)
    factor = sqrt(4*pi) .* (((2 .*ls) .+ 1).^(-0.5))
    Plm_arr = factor .* sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, cos(θ))[1:length(ls)]

    Ylm_vec = Vector{Complex{Float64}}(undef, (l_max+1)^2)
    Ylm_vec[1] = Plm_arr[1]

    ind1 = 1
    ind2 = 1
    for i = 1:l_max
        inds1 = (ind1+i):(ind1+2i)
        Ylm_vec[(ind2+i):(ind2+3i)] =
        [reverse(conj(Plm_arr[inds1[2:end]])); Plm_arr[inds1]]
        ind1 += i
        ind2 += 2i
    end

    ls, ms = spherical_harmonic_indices(l_max)
    Ylm_vec = exp.(ms .* (im*φ)) .* Ylm_vec

    return Ylm_vec
end

function anms(n::Int64,m::Int64)
    """
    The normalization factors for constructing moments.
    """
    den = sqrt(factorial(big(n-m))*factorial(big(n+m)))

    return (-1)^n/den
end

function Inm(n::Int64,m::Int64,p::Int64,ρ::Float64,θ::Float64,ϕ::Float64)
    """
    The inner functions for constructing moments.
    """
    Ynm = spherical_harmonics(p,θ,ϕ)
    ind = lm_to_spherical_harmonic_index(n,m)

    return im^float(-abs(m))*anms(n,m)*Ynm[ind]*ρ^n

end

function Onm(n::Int64,m::Int64,p::Int64,ρ::Float64,θ::Float64,ϕ::Float64)
    """
    The outer functions for constructing moments.
    """
    Ynm = spherical_harmonics(p,θ,ϕ)
    ind = lm_to_spherical_harmonic_index(n,m)
    num = ((-1)^n)*im^float(abs(m))*Ynm[ind]
    den = anms(n,m)*ρ^(n+1)

    return num/den

end

end #module
