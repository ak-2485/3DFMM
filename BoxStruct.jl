"""
Box structures for the 3D FMM and FMM functions on boxes.
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
    level::Int64
    index::Int64
    parent::Int64
    children::Set{Int64}
    particles::Set{Int64}
    numparticles::Int64
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
    returns true if box2 is a colleague of box1. Colleagues are on the same
    level and share a vertex.
    """
    A = Set(box1.vertices)
    B = Set(box2.vertices)
    return box1.level == box2.level && !isempty(intersect(A,B))
end

function arewellseparated(box1::Box, box2::Box)
    """
    Returns true if box1 and box2 are well separated: they are on the same level
    and aren't colleagues (i.e. they don't share a vertex).
    """
    return box1.level == box2.level && !iscolleague(box1,box2)
end

function areadjacent(box1::Box, box2::Box)
    """
    Returns true if box2 is adjacent to box1
    """
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
        Returns the set of integer ids of the particles in the grid dictionary of particles
        that are in the bounds of box1.
        """
    particleset = Set{Int}()
    minBound = box1.min_bound
    maxBound = box1.max_bound
    for partid in keys(particles)
        (x,y,z) = particles[partid][1]
        if x < maxBound[1] && y < maxBound[2] && z < maxBound[3] &&
            x >= minBound[1] && y >= minBound[2] && z >= minBound[3]
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
    ρxyz < 1e-6 ? θ = pi/2 : θ = acos(z/ρxyz)

    return ρxyz, θ, ϕ
end

function mcoeftrans(box1::Box, box2::Box, p::Int64)
    """
    Returns a vector of multipole coefficients translated from box1 to box2.
    """

    Mnm = box1.multipole_coef
    ρ,α,β = spherecenter(box2, box1.center)
    Mjk = zeros(ComplexF64,length(Mnm))

    for j = 0:p
        for k = -j:j
            ind1 = lm_to_spherical_harmonic_index(j,k)
            for n = 0:j
                for m = -n:n
                    if abs(k-m) <= j-n
                        ind2 = lm_to_spherical_harmonic_index(n,m)
                        Mjk[ind1] += Mnm[ind2] * Inm(j-n,k-m,p,ρ,α,β)
                    end
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
    Mnm = box1.multipole_coef
    ρ,α,β = spherecenter(box2, box1.center)
    Ljk = zeros(ComplexF64,length(Mnm))

    for j = 0:p
        for k = -j:j
            ind = lm_to_spherical_harmonic_index(j,k)
            for n = 0:p
                for m = -n:n
                    ind1 = lm_to_spherical_harmonic_index(n,m)
                    Ljk[ind] += Mnm[ind1] * Onm(j+n,-k-m,2p,ρ,α,β)
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
    Lnm = box1.local_coef
    ρ,α,β = spherecenter(box2, box1.center)
    Ynm = spherical_harmonics(p,α,β)
    Ljk = zeros(ComplexF64,length(Lnm))

    for j = 0:p
        for k = -j:j
            ind = lm_to_spherical_harmonic_index(j,k)
            for n = j:p
                for m = -n:n
                    if abs(m-k) <= n-j
                        ind2 = lm_to_spherical_harmonic_index(n,m)
                        Ljk[ind] += Lnm[ind2] * Inm(n-j,m-k,p,ρ,α,β)
                    end
                end#ms
            end#ns
        end#ks
    end#js

    return Ljk
end #lcoeftrans


end #end module
