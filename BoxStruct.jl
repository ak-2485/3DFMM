
"""
"""
module BoxStruct

using CoefHelpers
using GSL

export Box, iscolleague, arewellseparated, areadjacent, particlesin, haschildren
export mcoeftrans, mtolconversion, lcoeftrans

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
    multipole_coef::Vector{Complex{Float64}}
    local_coef::Vector{Complex{Float64}}
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
        this.multipole_coef = Vector{Complex{Float64}}()
        this.local_coef = Vector{Complex{Float64}}()
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
        if x <= maxBound[1] && y < maxBound[2] && z <= maxBound[3] &&
            x > minBound[1] && y >= minBound[2] && z > minBound[3]
            push!(particleset,partid)
        end
    end
    return particleset
end

function haschildren(box1::Box)
    """
    Returns true if box1 has children.
    """
    return !isempty(box1.children)
end

function spherecenter(box::Box, origin::Array{Float64,1})
    """
    Yields the spherical coordinates of the center of "box" with respect to the
    point "origin"
    """
    x = box.center[1]
    y = box.center[2]
    z = box.center[3]
    cx = origin[1]
    cy = origin[2]
    cz = origin[3]
    (x,y,z) = (x-cx,y-cy,z-cz)
    ρxy = sqrt(x^2 + y^2)
    ρxyz = sqrt(x^2 + y^2 + z^2)
    ϕ = atan(y,x)
    #ρxy < 1e-6 ? ϕ = 0.0 : ϕ = acos(x/ρxy)
    ρxyz < 1e-6 ? θ = pi/2 : θ = acos(z/ρxyz)
    #θ = atan(ρxy,z)

    return ρxyz, θ, ϕ
end

function anms(n::Int64,m::Int64)
    """
    """
    den = sqrt(factorial(big(n-m))*factorial(big(n+m)))

    return 1/den
end

function mcoeftrans(box1::Box, box2::Box, p::Int64)
    """
    Returns a vector of multipole coefficients translated from box1 to box2.
    """

    println("child center: ", box1.center)
    println("parent center: ", box2.center)

    ns,ms = spherical_harmonic_indices(p)
    len = length(ns)
    Onm = box1.multipole_coef
    ρ,α,β = spherecenter(box1, box2.center)
    Ynm = spherical_harmonics(p,α,β)
    Mjk = zeros(ComplexF64,len)

    for j = 0:p
        for k = -j:j
            ind = lm_to_spherical_harmonic_index(j,k)
            for n = 0:j
                for m = max(k+n-j,-n):min(k+j-n,n)
                    ind2 = lm_to_spherical_harmonic_index(j-n,k-m)
                    ind3 = lm_to_spherical_harmonic_index(n,-m)
                    num = Onm[ind2] * (-1)^m * anms(n,m) * anms(j-n,k-m) * ρ^n * Ynm[ind3]
                    Mjk[ind] += num/anms(j,k)
                end#ms
            end#ns
        end#ks
    end#js

    return Mjk
end #mcoeftrans

function mtolconversion(box1::Box, box2::Box, p::Int64)
    """
    Return a vector of the multipole coefficients of box1 converted to local
    coefficients centered at box2.
    """
    Onm = box1.multipole_coef
    ρ,α,β = spherecenter(box1, box2.center)
    Ynm = spherical_harmonics(2p,α,β)

    Ljk = zeros(ComplexF64,(p+1)^2)

    for j = 0:p
        for k = -j:j
            ind = lm_to_spherical_harmonic_index(j,k)
            for n = 0:p
                for m = -n:n
                    ind1 = lm_to_spherical_harmonic_index(n,m)
                    ind2 = lm_to_spherical_harmonic_index(j+n,m-k)
                    num = Onm[ind1] * (-1)^(n+k) * anms(n,m) * anms(j,k) * Ynm[ind2]
                    den = anms(j+n,m-k) * ρ^(j+n+1)
                    Ljk[ind] += num/den
                end#ms
            end#ns
        end#ks
    end#js

    return Ljk
end


function lcoeftrans(box1::Box, box2::Box, p::Int64)
    """
    Return a vector of the local coefficients for box1 translated to the center
    of box2.
    """

    ns,ms = spherical_harmonic_indices(p)
    Onm = box1.local_coef
    ρ,α,β = spherecenter(box1, box2.center)
    Ynm = spherical_harmonics(p,α,β)
    len = length(ns)
    Ljk = zeros(ComplexF64,len)

    for j = 0:p
        for k = -j:j
            ind = lm_to_spherical_harmonic_index(j,k)
            for n = j:p
                for m = k-n+j:k-j+n
                    ind1 = lm_to_spherical_harmonic_index(n-j,m-k)
                    ind2 = lm_to_spherical_harmonic_index(n,m)
                    num = Onm[ind2] * (-1)^(n-j) * anms(j,k) * anms(n-j,m-k) * ρ^(n-j) * Ynm[ind1]
                    den = anms(n,m)
                    Ljk[ind] += num/den
                end#ms
            end#ns
        end#ks
    end#js

    return Ljk
end #lcoeftrans


end #end module
