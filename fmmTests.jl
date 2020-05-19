


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


function fmmTest()

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
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
nmax = 1000
p = 5
println("nmax:", nmax)
println("p:", p)
num = 2000
println("numparticles:", num)

minbound = [-0.5,-0.5,-0.5]
maxbound = [0.5,0.5,0.5]

coords = rand(Uniform(-0.5,0.5), num,3)
sz=size(coords)[1]

particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()
for i = 1:sz
    x,y,z = coords[i,:]
    coord = (x,y,z)
    q = rand(-5:0.1:5)
    particledict[i] = (coord,q,0.0)
    particles[i] = (coord,q)
end

@time begin
ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)
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

randomTest1()

function uniformgridTest()
println("Testing function fast mulipole: Uniform grid")
nmax = 1000
p = 10
println("maximum number of particles in a leaf box:", nmax)
println("expansion order:", p)

coordsx = [x for x in 0.25:0.025:0.75 for y in 0.25:0.025:0.75 for z in 0.25:0.025:0.75]
coordsy = [y for x in 0.25:0.025:0.75 for y in 0.25:0.025:0.75 for z in 0.25:0.025:0.75]
coordsz = [z for x in 0.25:0.025:0.75 for y in 0.25:0.025:0.75 for z in 0.25:0.025:0.75]
sz=length(coordsx)
coords = zeros(sz,3)
coords[:,1] = coordsx
coords[:,2] = coordsy
coords[:,3] = coordsz

minbound = [0.0,0.0,0.0]
maxbound = [1.0,1.0,1.0]

particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()
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
ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)
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
#    println(vals[3])
end
#println(particledict)
println(E)
end #test2

uniformgridTest()
