
"""
"""
module GridStruct

using BoxStruct
using CoefHelpers
using LinearAlgebra

export Grid, carttosphere, mcoef!, lcoef!, sphcoords, carttospherepoint
export localatpoints!, directatpoints!, multipoleatpoints!

mutable struct Grid
    """
    Represents a box in the 3D FMM.

    fields:
    """
    boxes::Dict{Int64,Box}
    particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}
    numparticles::Int64
    levels::Dict{Int64,Set{Int64}}
    numlevels::Int64
    minbound::Array{Float64,1}
    maxbound::Array{Float64,1}
    center::Array{Float64,1}
    size::Int64
    signarray::Array{Array{Int64,1},1}
    nmax::Int64
    function Grid()
        this = new()
        this.boxes = Dict{Int64,Box}()
        this.particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}()
        this.numparticles = 0
        this.levels = Dict{Int64,Set{Int64}}()
        this.numlevels = 0
        this.minbound = zeros(Float64,3)
        this.maxbound = zeros(Float64,3)
        this.center = zeros(Float64,3)
        this.size = 1
        # An array for calculating vertices and centers of refinement boxes
        # arranged s.t. |7,8|5,6|
        #              |1,2|3,4|
        this.signarray = [[-1,-1,1],[-1,-1,-1],[1,-1,1],[1,-1,-1],
                    [1,1,1],[1,1,-1],[-1,1,1],[-1,1,-1]]
        this.nmax = 2
        @assert this.nmax > 1 ["minimum number of particles in box must be greater than 1"]
        return this
    end
end

function sphcoords(particleid::Int64, grid::Grid)
    """
    """

    particle = grid.particles[particleid]
    (x,y,z) = particle[1]
    ρxy = sqrt(x^2 + y^2)
    ρxyz = sqrt(x^2 + y^2 + z^2)
    ϕ = atan(y,x)
    #ρxy < 1e-6 ? ϕ = 0.0 : ϕ = acos(x/ρxy)
    ρxyz < 1e-6 ? θ = pi/2 : θ = acos(z/ρxyz)
    #θ = atan(ρxy,z)

    return ρxyz, θ, ϕ

end

function carttosphere(particleid::Int64, box1::Box, grid::Grid)
    """
    Yields the spherical coordinates for a particle
    with respect to the center of box1.

    inputs:

    preconditions: particleid is a key in particles.
    returns: a triple of the spherical coordinates of the particles with
        respect to the center of box1.
    """
    # the center of the box
    cx = box1.center[1]
    cy = box1.center[2]
    cz = box1.center[3]
    particle = grid.particles[particleid]
    (x,y,z) = particle[1]
    (x,y,z) = (x-cx,y-cy,z-cz)
    ρxy = sqrt(x^2 + y^2)
    ρxyz = sqrt(x^2 + y^2 + z^2)
    ϕ = atan(y,x)
    #ρxy < 1e-6 ? ϕ = 0.0 : ϕ = acos(x/ρxy)
    ρxyz < 1e-6 ? θ = pi/2 : θ = acos(z/ρxyz)
    #θ = atan(ρxy,z)

    return ρxyz, θ, ϕ
end

function carttospherepoint(point::Array{Float64,1}, origin::Array{Float64,1}, grid::Grid)
    """
    Yields the spherical coordinates for a point
    with respect to the center of box1.

    inputs:

    preconditions: particleid is a key in particles.
    returns: a triple of the spherical coordinates of the particles with
        respect to the center of box1.
    """
    # the center of the box
    cx = origin[1]
    cy = origin[2]
    cz = origin[3]
    (x,y,z) = point
    (x,y,z) = (x-cx,y-cy,z-cz)
    ρxy = sqrt(x^2 + y^2)
    ρxyz = sqrt(x^2 + y^2 + z^2)
    ϕ = atan(y,x)
    #ρxy < 1e-6 ? ϕ = 0.0 : ϕ = acos(x/ρxy)
    ρxyz < 1e-6 ? θ = pi/2 : θ = acos(z/ρxyz)
    #θ = atan(ρxy,z)

    return ρxyz, θ, ϕ
end

function mcoef!(box::Box, p::Int64, grid::Grid)
    """
    Find the multipole expansion coefficients Mnm for a given box, about its center, where
        M_n^m =
        sum_i=1^numparticles [q_i * r_i^n * Y_n^-m(theta_i,phi_i)]

    inputs:
    returns: A vector of the coefficints M_nm for a given box. E.g.
        [M_0^0 M_1^-1 M_1^0  M_1^1... M_p^p]
    """
    particles = box.particles
    Mnm = zeros(ComplexF64,length(box.multipole_coef))
    for n = 0:p
        for m = -n:n
            for partid in particles
                particle = grid.particles[partid]
                q = particle[2]
                ρ, θ, ϕ = carttosphere(partid,box,grid)
                Ynm = spherical_harmonics(p, θ, ϕ)
                ind1 = lm_to_spherical_harmonic_index(n,m)
                ind2 = lm_to_spherical_harmonic_index(n,-m)
                Mnm[ind1] += (-1)^m *(q*ρ^n)*Ynm[ind2]
            end
        end
    end

    box.multipole_coef .+= Mnm

end # fun

function lcoef!(box1::Box, box2::Box, p::Int64, grid::Grid)
    """
    Determine the contributions to box1's local expansion from particles in
    box 2.
    """
    Ljk = zeros(ComplexF64,(p+1)^2)
    for j = 0:p
        for k = -j:j
            for partid in box2.particles
                particle = grid.particles[partid]
                q = particle[2]
                ρ, θ, ϕ = carttosphere(partid,box2,grid)
                Yjk = spherical_harmonics(p, θ, ϕ)
                ind1 = lm_to_spherical_harmonic_index(j,k)
                ind2 = lm_to_spherical_harmonic_index(j,-k)
                Ljk[ind1] += (q/ρ^(j+1))*Yjk[ind2]
            end
        end
    end

    box1.local_coef = Ljk

end#function

function localatpoints!(box::Box, grid::Grid, p::Int64)
    """
    Updates the potential for all particles in box from box's local expansion.
    """

    Ljk = box.local_coef
    for particleid in box.particles
        ψ = 0
        ρ, θ, ϕ = carttosphere(particleid, box, grid)
        Yjk = spherical_harmonics(p, θ, ϕ)
        for j = 0:p
            for k = -j:j
                ind1 = lm_to_spherical_harmonic_index(j,k)
                ψ += Ljk[ind1] * Yjk[ind1] * ρ^j
            end
        end
        coords,q,pot = grid.particles[particleid]
        vals = (coords, q, pot + ψ)
        grid.particles[particleid] = vals

    end

end #evalpotential

function directatpoints!(box1::Box, box2::Box, grid::Grid, p::Int64)
    """"
    Updates the potential for all particles in box1 by including directly
    calculated contributions from box2.
    """

    for particleid1 in box1.particles
        for particleid2 in box2.particles
            if particleid1 != particleid2
                coords,q,ϕ = grid.particles[particleid1]
                coordsi = grid.particles[particleid2][1]
                qi = grid.particles[particleid2][2]
                den = norm(coords.-coordsi)
                vals = (coords, q, ϕ + qi/den)
                grid.particles[particleid1] = vals
            end
        end
    end

end #function

function multipoleatpoints!(box1::Box, box2::Box, grid::Grid, p::Int64)
    """
    Updates the potential for all particles in box1 by including contributions
    from box2's multipole expansion.
    """

    Mnm = box2.multipole_coef
    for particleid in box1.particles
        ψ = 0
        ρ, θ, ϕ = carttosphere(particleid, box1, grid)
        Ynm = spherical_harmonics(p, θ, ϕ)
        for n = 0:p
            for m = -n:n
                ind1 = lm_to_spherical_harmonic_index(m,n)
                ψ += Mnm[ind1] * Ynm[ind1] * ρ^-(n+1)
            end
        end
        coords,q,pot = grid.particles[particleid]
        vals = (coords, q, pot + ψ)
        grid.particles[particleid] = vals
    end

end #evalpotential

end  # module GridStruct
