
"""
"""
module BoxStruct

using CoefHelpers
using GSL

export Box, iscolleague, arewellseparated, areadjacent, particlesin, haschildren
export mcoeftrans, mtolconversion, lcoefsum!, lcoeftrans

mutable struct Box
    """
    Represents a box in the 3D FMM.

    fields:
    """
    # Boxes have 20 fields
    center::Array{Float64,1}
    size::Float64
    vertices::Array{Array{Float64,1},2}
    min_bound::Array{Float64,1}
    max_bound::Array{Float64,1}
    level::Int
    index::Int
    parent::Int
    children::Set{Int64}
    particles::Set{Int64}
    numparticles::Int
    colleagues::Set{Int64}
    L1::Set{Int64}
    L2::Set{Int64}
    L3::Set{Int64}
    L4::Set{Int64}
    uplist::Set{Int64}
    downlist::Set{Int64}
    northlist::Set{Int64}
    southlist::Set{Int64}
    eastlist::Set{Int64}
    westlist::Set{Int64}
    multipole_coef::Array{ComplexF64,1}
    local_coef::Array{ComplexF64,1}
    function Box()
        this = new()
        this.center = zeros(Float64,1)
        this.vertices = fill(zeros(Float64,3),1,8)
        this.max_bound = zeros(Float64,3)
        this.min_bound = zeros(Float64,3)
        this.size = 1
        this.level =0
        this.parent = 0
        this.children=Set{Int64}()
        this.index= 0
        this.particles=Set{Int64}()
        this.numparticles=0
        this.colleagues = Set{Int64}()
        this.L1 = Set{Int64}()
        this.L2 = Set{Int64}()
        this.L3 = Set{Int64}()
        this.L4 = Set{Int64}()
        this.uplist = Set{Int64}()
        this.downlist = Set{Int64}()
        this.northlist = Set{Int64}()
        this.southlist = Set{Int64}()
        this.eastlist = Set{Int64}()
        this.westlist = Set{Int64}()
        this.multipole_coef = Array{ComplexF64,1}()
        this.local_coef = Array{ComplexF64,1}()
        return this
    end
end

function iscolleague(box1::Box, box2::Box)
    """
    """
    # Colleagues are on the same level and share a vertex
    A = Set(box1.vertices)
    B = Set(box2.vertices)
    return box1.level == box2.level && !isempty(intersect(A,B))
end

function arewellseparated(box1::Box, box2::Box)
    """
    """
    # Well separated boxes are on the same level and aren't colleagues
    # i.e. they don't share a vertex
    return box1.level == box2.level && !iscolleague(box1,box2)
end

function areadjacent(box1::Box, box2::Box)
    """
    """
    # https://stackoverflow.com/questions/5009526/overlapping-cubes
    A = Set(box1.vertices)
    B = Set(box2.vertices)
    if !isempty(intersect(A,B)) return true end
    box1Max = box1.max_bound
    box1Min = box1.min_bound
    box2Max = box2.max_bound
    box2Min = box2.min_bound
    checkX = box1Max[1] < box2Min[1] || box2Max[1] < box1Min[1]
    checkY = box1Max[2] < box2Min[2] || box2Max[2] < box1Min[2]
    checkZ = box1Max[3] < box2Min[3] || box2Max[3] < box1Min[3]
    return !(checkX || checkY || checkZ)
end

function particlesin(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}},
        box1::Box)
        """
        """
    particleset = Set{Int}()
    minBound = box1.min_bound
    maxBound = box1.max_bound
    for partid in keys(particles)
        (x,y,z) = particles[partid][1]
        if x < maxBound[1] && y < maxBound[2] && z < maxBound[3] &&
            x > minBound[1] && y > minBound[2] && z > minBound[3]
            push!(particleset,partid)
        end
    end
    return particleset
end

function haschildren(box1::Box)
    """
    """
    return !isempty(box1.children)
end

function spherecenter(box::Box)
    """
    Yields the spherical coordinates of the center of a box
    """
    cx = box.center[1]
    cy = box.center[2]
    cz = box.center[3]
    ρ = sqrt((cx)^2 + (cy)^2 + (cz)^2)
    θ = atan((cy)/(cx))
    ϕ = acos((cz)/ρ)

    return ρ, θ, ϕ
end

function mcoeftrans(box1::Box, p::Int64)
    """
    Return a vector of translated multipole coefficients of box1.
    """

    ns,ms = spherical_harmonic_indices(p)
    len = length(ns)

    Anm = zeros(ComplexF64,len)
    for i = 1:len
        ind = lm_to_spherical_harmonic_index(ns[i],ms[i])
        n = ns[i]
        m = ms[i]
        den = sqrt(factorial(n-m)*factorial(n+m))
        Anm[ind] = 1/den
    end

    Onm = box1.multipole_coef
    ρ,α,β = spherecenter(box1)
    Ynm = spherical_harmonics(p,α,β)
    #ρn = (ρ .* ones(len)).^ns

    Mjk = zeros(ComplexF64,len)
    Ys = Vector{ComplexF64}()
    As1 = Vector{ComplexF64}()
    As2 = Vector{ComplexF64}()
    Js = Vector{ComplexF64}()
    Os = Vector{ComplexF64}()
    ρn = Vector{ComplexF64}()

    for j = 0:p
        for k = -j:j
            for n = 0:j
                for m = max(k+n-j,-n):min(k+j-n,n)
                    ρn = [ρn; ρ^n]
                    ind1 = lm_to_spherical_harmonic_index(n,-m)
                    Ys = [Ys ; Ynm[ind1]]
                    ind2 = lm_to_spherical_harmonic_index(j-n,k-m)
                    As1 = [As1 ; Anm[ind2]]
                    Os = [Os; Onm[ind2]]
                    ind3 = lm_to_spherical_harmonic_index(n,m)
                    As2 = [As2 ; Anm[ind3]]
                    jterm = (-1)^m
                    Js = [Js ; jterm]
                end#ms
            end#ns
            ind = lm_to_spherical_harmonic_index(j,k)
            terms = Js .* As2 .* As1 .* ρn .* Ys
            Mjk[ind] = (1/Anm[ind])*(Os' * terms)
        end#ks
    end#js

    return Mjk
end #mcoeftrans

function mtolconversion(box1::Box, p::Int64)
    """
    Return a vector of the multipole coefficients for box1 converted to local
    coefficients.
    """
    ns,ms = spherical_harmonic_indices(2p)
    len = length(ns)
    sz = (p+1)^2

    Anm = zeros(ComplexF64,len)
    for i = 1:len
        ind = lm_to_spherical_harmonic_index(ns[i],ms[i])
        n = ns[i]
        m = ms[i]
        den = sqrt(factorial(n-m)*factorial(n+m))
        Anm[ind] =1/den
    end

    Onm = box1.local_coef
    ρ,α,β = spherecenter(box1)
    Ynm = spherical_harmonics(2p,α,β)

    Ljk = zeros(ComplexF64,sz)
    Ys = Vector{ComplexF64}(undef,sz)
    As1 = Vector{ComplexF64}(undef,sz)
    As2 = Vector{ComplexF64}(undef,sz)
    Js = Vector{ComplexF64}(undef,sz)
    Os = Vector{ComplexF64}(undef,sz)
    ρn = Vector{ComplexF64}(undef,sz)

    for j = 0:p
        for k = -j:j
            for n = 0:p
                for m = -n:n
                    ind1 = lm_to_spherical_harmonic_index(n,m)
                    ρn[ind1] = ρ^-(j+n+1)
                    Js[ind1] = (-1)^(n+k)
                    ind2 = lm_to_spherical_harmonic_index(j+n,m-k)
                    Ys[ind1] = Ynm[ind2]
                    As1[ind1] = 1/Anm[ind2]
                end#ms
            end#ns
            ind = lm_to_spherical_harmonic_index(j,k)
            terms = Anm[1:sz] .* Anm[ind] .* Js .* As1 .* ρn .* Ys
            Ljk[ind] = (Onm' * terms)
        end#ks
    end#js

    return Ljk
end

function lcoefsum!(box1::Box, coefficients::Vector{ComplexF64}, p::Int64)
    """
    Add local expansion coefficients, "coefficients", to the local coefficients of box1.
    """
    # Add coefficients from coefvec to box1's local expansion coefficients
    box1.local_coef .+= coefficients
end

function lcoeftrans(box1::Box, p::Int64)
    """
    Return a vector of the translated local coefficients for box1.
    """
    ns,ms = spherical_harmonic_indices(p)
    len = length(ns)

    Anm = zeros(ComplexF64,len)
    for i = 1:len
        ind = lm_to_spherical_harmonic_index(ns[i],ms[i])
        n = ns[i]
        m = ms[i]
        den = sqrt(factorial(n-m)*factorial(n+m))
        Anm[ind] = 1/den
    end

    Onm = box1.local_coef
    ρ,α,β = spherecenter(box1)
    Ynm = spherical_harmonics(p,α,β)
    #ρn = (ρ .* ones(len)).^ns

    Ljk = zeros(ComplexF64,len)
    Ys = Vector{ComplexF64}()
    As1 = Vector{ComplexF64}()
    As2 = Vector{ComplexF64}()
    Js = Vector{ComplexF64}()
    Os = Vector{ComplexF64}()
    ρn = Vector{ComplexF64}()

    for j = 0:p
        for k = -j:j
            for n = j:p
                for m = k-n+j:k-j+n
                    ρn = [ρn; ρ^(n-j)]
                    ind1 = lm_to_spherical_harmonic_index(n-j,m-k)
                    Ys = [Ys ; Ynm[ind1]]
                    As1 = [As1 ; Anm[ind1]]
                    ind2 = lm_to_spherical_harmonic_index(n,m)
                    Os = [Os; Onm[ind2]]
                    As2 = [As2 ; 1/Anm[ind2]]
                    jterm = (-1)^(n-j)
                    Js = [Js ; jterm]
                end#ms
            end#ns
            ind = lm_to_spherical_harmonic_index(j,k)
            terms = Js .* As2 .* As1 .* ρn .* Ys
            Ljk[ind] = Anm[ind]*(Os' * terms)
        end#ks
    end#js

    return Ljk
end #lcoeftrans


end #end module
