
module CoefHelpers

using GSL

export spherical_harmonics, spherical_harmonic_indices
export lm_to_spherical_harmonic_index, anms, Inm, Onm
"""
"""
# https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/src/specialfunctions.jl
function spherical_harmonic_indices(lmax::Int)
    ls = [l for l in 0:lmax for m in -l:l]
    ms = [m for l in 0:lmax for m in -l:l]

    return ls, ms
end

function associated_legendre_indices(lmax::Int)
    ls = [l for l in 0:lmax for m in 0:l]
    ms = [m for l in 0:lmax for m in 0:l]

    return ls, ms
end

lm_to_spherical_harmonic_index(l::Int,m::Int)::Int = l^2 + m + l + 1

function spherical_harmonics(l_max::Int, θ::T, φ::T) where T <: AbstractFloat
    """
    `spherical_harmonics(l_max::Int, θ::T, φ::T)`
    returns a vector of all spherical harmonics with degree `l <= l_max`.
    The degree and order (indices) of the elements of the vector are given by
    `spherical_harmonics_indices(l_max::Int)`.
    The associated legendre polynomials are taken from the package GSL.jl.
    """

    ls, ms = associated_legendre_indices(l_max)
    factor = sqrt(4*pi) .* (((2 .*ls) .+ 1).^(-0.5))
    Plm_arr = factor .* sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, cos(θ))[1:length(ls)]

    Ylm_vec = Vector{Complex{T}}(undef, (l_max+1)^2)
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
    """
    den = sqrt(factorial(big(n-m))*factorial(big(n+m)))

    return (-1)^n/den
end

function Inm(n::Int64,m::Int64,p::Int64,ρ::Float64,θ::Float64,ϕ::Float64)
    """
    """
    Ynm = spherical_harmonics(p,θ,ϕ)
    ind = lm_to_spherical_harmonic_index(n,m)

    return im^float(-abs(m))*anms(n,m)*Ynm[ind]*ρ^n

end

function Onm(n::Int64,m::Int64,p::Int64,ρ::Float64,θ::Float64,ϕ::Float64)
    """
    """
    Ynm = spherical_harmonics(p,θ,ϕ)
    ind = lm_to_spherical_harmonic_index(n,m)
    num = ((-1)^n)*im^float(abs(m))*Ynm[ind]
    den = anms(n,m)*ρ^(n+1)

    return num/den

end

function wigner_matrices_indices(lmax::Int)
    ls = [l for l in 0:lmax for m in -l:l for k in -l:l]
    ms = [m for l in 0:lmax for m in -l:l for k in -l:l]
    ks = [k for l in 0:lmax for m in -l:l for k in -l:l]

    return ls, ms, ks
end

# Order is [l,k,m]
lmk_to_wigner_index(l::Int,k::Int,m::Int)::Int = ((((4l + 6)l + 6m + 5)l + 3(m +k)) / 3)+1

function sgn(ms::Array{Int64,1})
    for i=1:length(ms)
        ms[i] < 0 ? ms[i]=-1 : ms[i]= 1
    end
    return ms
end

# https://arxiv.org/pdf/1403.7698.pdf
function wigner_matrices(lmax::Int, θ::T) where T <: AbstractFloat
    """
    Following eq. 41
    """
    # Calculate the coefficients
    ls, ms = spherical_harmonic_indices(2lmax)
    len = length(ls)
    inds1 = 1:len

    blm_vec = Vector{Complex{Float64}}(undef, len)
    num = (ls[inds1] .- ms[inds1] .- 1)*(ls[inds1] .- ms[inds1])
    den = (2ls[inds1] .- 1)*(2ls[inds1] .+ 1)
    bsqrtlsms = (num./den).^0.5
    blm_vec[inds1] = sgn(ms[inds1]) .* bsqrtlsms

    alm_vec = Vector{Complex{Float64}}(undef, len)
    for i = inds1
        l = ls[i]
        m = ms[i]
        j = lm_to_spherical_harmonic_index(l,m)
        if l < abs(m)
            alm_vec[j] = 0
        else
            num = (l + 1 + m) * (l + 1 - m)
            den = (2l + 1) * (2l + 3)
            alm_vec[j] = (num/den)^0.5
        end
    end

    Ylm_theta_vec = spherical_harmonics_theta(2lmax, θ)

    ls, ks, ms = wigner_matrices_indices(2lmax)
    len = length(ls)

    Hlkm_vec = Vector{Complex{Float64}}(undef, len)
    for i = 1:len
        if ms[i] == 0
            ind1 = lkm_to_wigner_index(ls[i],ks[i],ms[i])
            ind2 = lkm_to_wigner_index(ls[i],ms[i],ks[i])
            ind3 = lm_to_spherical_harmonic_index(ls[i],ks[i])
            Dlkm[ind1] = Ylm_theta_vec[ind3]
            # D_l^{k}{0} = D_l^{0}{k}
            Dlkm[ind2] = Dlkm[ind1]
        end
    end
    for n = 2lmax:-1:2
    lst = min(2lmax-n,n-2)
        for m = 0:lst
            for k = -n+1:n-1
                ind1 = lkm_to_wigner_index(n-1,k,m+1)
                indxb1 = lm_to_spherical_harmonic_index(n,m)
                invb = 1/blm_vec[indxb1]

                ind2 = lkm_to_wigner_index(n,k+1,m)
                indxb2 = lm_to_spherical_harmonic_index(n,-k-1)
                term2 = 0.5 * (1 - cos(θ))*blm_vec[indxb2]*Dlkm[ind2]

                ind3 = lkm_to_wigner_index(n,k-1,m)
                indxb3 = lm_to_spherical_harmonic_index(n,k-1)
                term3 = 0.5 * (1 + cos(θ))*blm_vec[indxb3]*Dlkm[ind3]

                ind4 = lkm_to_wigner_index(n,k,m)
                indxa4 = lm_to_spherical_harmonic_index(n-1,k)
                term4 = -sin(θ)*alm_vec[indxa4]*Dlkm[ind4]

                Dlkm[ind1] = term2 - term3 - term4
                Dlkm[ind1] = invb * Dlkm[ind1]

                # apply symmetry
                ind5 = lkm_to_wigner_index(n-1,m+1,k)
                mind1 = lkm_to_wigner_index(n-1,-k,-(m+1))
                mind5 = lkm_to_wigner_index(n-1,-(m+1),-k)
                Dlkm[ind5] = Dlkm[ind1]
                Dlkm[mind1] = Dlkm[ind1]
                Dlkm[mind5] = Dlkm[ind5]
            end
        end
    end
end #function wigner_matrices



end #module
