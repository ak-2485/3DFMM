
using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using Random, Distributions
import FMMSteps
import FMM3D
import BoxStruct
import GridStruct

particles = Dict([(1,((0.24,0.22,0.26),3.0)),(2,((0.72,0.76,0.72),4.0)),
(3,((0.54,0.25,0.25),3.0)),(4,((0.73,0.65,0.72),2.0)), (5,((0.12,0.22,0.22),2.0)),
(6,((0.72,0.72,0.14),3.0))])


function mcoefTest()
    """
    Test base case multipole expansions: for box 0, entire domain.
    """
    println("Testing multipole expansion")
    p = 2
    nmin = 8 # just has to be greater than particle list tested to get single box
    grid = FMM3D.step0(particles,nmin,p)
    box = grid.boxes[0]

    Mnm = zeros(ComplexF64,(p+1)^2) #initialize; contains zeros at first

    # Calculate multipoles directly
    for n = 0:p
        for m = -n:n
            for (partid, vals) in particles
                q = vals[2]
                ρ, θ, ϕ = sphcoords(partid,grid)
                Ynm = spherical_harmonics(p, θ, ϕ)
                ind1 = lm_to_spherical_harmonic_index(n,m)
                ind2 = lm_to_spherical_harmonic_index(n,-m)
                Mnm[ind1] += (q*ρ^n)*Ynm[ind2]
            end
        end
    end


    # Cacluate multipole coefficients used in fastmultipole function
    FMM3D.step1(grid,p)
    num = 0
    den = 0
    for i = 1:length(Mnm)
        Mnmexpected = Mnm[i]
        Mnmfun = box.multipole_coef[i]
        num += abs(Mnmexpected - Mnmfun)^2
        den += abs(Mnmexpected)^2
    end
    E = num/den
    @test E < 10^-6
end

mcoefTest()

function mtranslationTest()
    p = 2
    nmin = 2
    grid = FMM3D.step0(particles,nmin,p)
    FMM3D.step1(grid,p)
    # get a leaf box
    child = grid.boxes[2]
    parent = grid.boxes[child.parent]
    # set its multipole coefficients
    GridStruct.mcoef!(child,p,grid)
    Onm = child.multipole_coef

    function expected(parent,child,p)
        ρ,θ,ϕ = BoxStruct.spherecenter(child,parent.center)
        Ynm = spherical_harmonics(p,θ,ϕ,true)
        Mjk = zeros(ComplexF64,4)
        for j = 0:1
            for k = -j:j
                ind = lm_to_spherical_harmonic_index(j,k)
                for n = 0:j
                    for m = max(k+n-j,-n):min(k+j-n,n)
                        ind2 = lm_to_spherical_harmonic_index(j-n,k-m)
                        ind3 = lm_to_spherical_harmonic_index(n,m)
                        A1 = 1/sqrt(factorial(n-m)*factorial(n+m))
                        A2 = 1/sqrt(factorial((j-n)-(k-m))*factorial((j-n)+(k-m)))
                        A3 = 1/sqrt(factorial(j-k)*factorial(j+k))
                        Mjk[ind] += (Onm[ind2] * (-1)^m * A1 * A2 * ρ^n * Ynm[ind3])/A3
                    end
                end
            end
        end
        return Mjk
    end
    e1 = expected(parent,child,p)

    function returned(parent,child,grid,p)
        Mjk = BoxStruct.mcoeftrans(child,parent,p)
        return Mjk[1:4]
    end
    e2 = returned(parent,child,grid,p)

    println("Testing multipole translation")
    num = 0
    den = 0
    for i in length(e1)
        num += abs(e1[i] - e2[i])^2
        den += abs(e1[i])^2
    end
    E = sqrt(num/den)
    @test E < 10^-6

end

#mtranslationTest()

function mtolconversionTest()
    p = 2
    nmin = 3
    grid = FMM3D.step0(particles,nmin,p)
    FMM3D.step1(grid,p)
    # get a leaf box
    child = grid.boxes[2]
    # set its multipole coefficients
    GridStruct.mcoef!(child,p,grid)
    Onm = L2box.multipole_coef
    L2box = grid.boxes[6]

    Ljk = BoxStruct.mtolconversion(L2box,child,p)

    Ljkcheck = zeros(ComplexF64,4)
    ρ,θ,ϕ = BoxStruct.spherecenter(L2.box,child.center)
    Ynm = spherical_harmonics(p,θ,ϕ)

end

function fmmTest()

    particledict = Dict{Int64,Tuple{Tuple{Float32,Float32,Float32},Float32,Float32}}()
    for (key,value) in particles
        coords = value[1]
        q = value[2]
        particledict[key] = (coords,q,0.0)
    end

    ϕ = FMM3D.fastmultipole(particles,2,1)

    for particleid1 in keys(particles)
        for particleid2 in keys(particles)
            if particleid1 != particleid2
                coords,q,ψ = particledict[particleid1]
                coordsi = particledict[particleid2][1]
                qi = particledict[particleid2][2]
                den = norm(coords.-coordsi)
                vals = (coords, q, ψ + qi/den)
                particledict[particleid1] = vals
            end
        end
    end

    num = 0
    den = 0
    for particleid in keys(particles)
        ϕ1 = ϕ[particleid][3]
        ϕ2 = particledict[particleid][3]
        num += abs(ϕ2 - ϕ1)^2
        den += abs(ϕ2)^2
    end
    E = sqrt(num/den)

    println(ϕ)
    println(particledict)
    println(E)
    println("Testing function fast mulipole")
    #@test E < 10^-6

end #fmmTest

#fmmTest()

function randomTest1()
println("Testing function fast mulipole: unform random distribution")
nmin = 75
p = 9
println("nmin:", nmin)
println("p:", p)
num = 400
println("numparticles:", num)

coords = rand(Uniform(-0.5,0.5), num,3)
sz=size(coords)[1]

particledict = Dict{Int64,Tuple{Tuple{Float32,Float32,Float32},Float32,Float32}}()
particles = Dict{Int64,Tuple{Tuple{Float32,Float32,Float32},Float32}}()
for i = 1:sz
    x,y,z = coords[i,:]
    coord = (x,y,z)
    q = rand(-5:0.1:5)
    particledict[i] = (coord,q,0.0)
    particles[i] = (coord,q)
end

@time begin
ϕ = FMM3D.fastmultipole(particles,nmin,p)
end #time

@time begin
for particleid1 in keys(particledict)
    for particleid2 in keys(particledict)
        if particleid1 != particleid2
            coords,q,ψ = particledict[particleid1]
            coordsi = particledict[particleid2][1]
            qi = particledict[particleid2][2]
            den = norm(coords.-coordsi)
            vals = (coords, q, ψ + qi/den)
            particledict[particleid1] = vals
        end
    end
end
end #time

num = 0
den = 0
for particleid in keys(particledict)
    ϕ1 = ϕ[particleid][3]
    ϕ2 = particledict[particleid][3]
    num += abs(ϕ2 - ϕ1)^2
    den += abs(ϕ2)^2
end
E = sqrt(num/den)
for (keys,vals) in ϕ
    #println(vals[3])
end
#println(particledict[1])
println(E)
end #randtest1

#randomTest1()

function uniformgridTest()
println("Testing function fast mulipole: Uniform grid")
nmax = 10
p = 4
println("maximum number of particles in a leaf box:", nmax)
println("expansion order:", p)

coordsx = [x for x in 0.25:0.25:0.75 for y in 0.25:0.25:0.75 for z in 0.25:0.25:0.75]
coordsy = [y for x in 0.25:0.25:0.75 for y in 0.25:0.25:0.75 for z in 0.25:0.25:0.75]
coordsz = [z for x in 0.25:0.25:0.75 for y in 0.25:0.25:0.75 for z in 0.25:0.25:0.75]
sz=length(coordsx)
coords = zeros(sz,3)
coords[:,1] = coordsx
coords[:,2] = coordsy
coords[:,3] = coordsz

particledict = Dict{Int64,Tuple{Tuple{Float32,Float32,Float32},Float32,Float32}}()
particles = Dict{Int64,Tuple{Tuple{Float32,Float32,Float32},Float32}}()
for i = 1:sz
    x,y,z = coords[i,:]
    coord = (x,y,z)
    q = (-1)^i*(x+y)/z
    particledict[i] = (coord,q,0.0)
    particles[i] = (coord,q)
end

num = sz
println("numparticles:", num)

@time begin
ϕ = FMM3D.fastmultipole(particles,nmax,p)
end #time

@time begin
for particleid1 in keys(particledict)
    for particleid2 in keys(particledict)
        if particleid1 != particleid2
            coords,q,ψ = particledict[particleid1]
            coordsi = particledict[particleid2][1]
            qi = particledict[particleid2][2]
            den = norm(coords.-coordsi)
            vals = (coords, q, ψ + qi/den)
            particledict[particleid1] = vals
        end
    end
end
end #time

num = 0
den = 0
for particleid in keys(particledict)
    ϕ1 = ϕ[particleid][3]
    ϕ2 = particledict[particleid][3]
    num += abs(ϕ2 - ϕ1)^2
    den += abs(ϕ2)^2
end
E = sqrt(num/den)
for (keys,vals) in ϕ
    println(vals[3])
end
println(particledict)
println(E)
end #test2

uniformgridTest()
