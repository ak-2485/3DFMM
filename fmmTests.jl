using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
#using Plots
using EnsembleTests
import FMMSteps
import FMM3D
import BoxStruct
import GridStruct


function randomTest1()
    println("Testing function fast mulipole: unform random distribution")
    nmax = 200
    p = 4
    num = 2000
    particles, particledict, minbound, maxbound = rndmparticledist3(num)

    @time begin
        ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)
    end #time
    println("maximum number of particles in a leaf box:", nmax)
    println("expansion order:", p)
    println("number of particles:", length(particles))

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
    println("Error is: ", E)

end #randtest1

#randomTest1()

function uniformgridTest()

    particles, particledict, minbound, maxbound = particledist2()

    println("Testing function fast mulipole: Uniform grid")
    nmax = 10000
    p = 6
    println("maximum number of particles in a leaf box:", nmax)
    println("expansion order:", p)
    println("number of particles:", length(particles))

    @time begin
        ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)
    end #time

    @time begin
        counter = 0
        while counter < 500
            for particleid1 in keys(particledict)
                counter +=1
                for particleid2 in keys(particledict)
                    counter +=1
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
    println("Error is: ", E)

end #test2

#uniformgridTest()

function uniformgridsmallTest()

    particles, particledict, minbound, maxbound = particledist1ex()

    println("Testing function fast mulipole: Uniform grid")
    nmax = 30
    p = 6
    println("maximum number of particles in a leaf box:", nmax)
    println("expansion order:", p)
    println("number of particles:", length(particles))

    @time begin
        ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)
    end #time

    @time begin
        counter = 0
        while counter < 500
            for particleid1 in keys(particledict)
                counter +=1
                for particleid2 in keys(particledict)
                    counter +=1
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
    println("Error is: ", E)

end #test2

uniformgridsmallTest()
