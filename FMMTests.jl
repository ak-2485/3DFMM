"""
Tests for the 3D FMM.
"""

module FMMTests

using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using EnsembleTests
using CPUTime

import FMMSteps
import FMM3D
import BoxStruct
import GridStruct

function uniformTest(numparts::Int64,nmax::Int64,p::Int64,d::Float64)

    particles, particledict, minbound, maxbound = rndmparticledist3(numparts,d)

    println("Testing function fast multipole: Uniform grid")

    @time @CPUtime ϕ = FMM3D.fastmultipole(particles,minbound,maxbound,nmax,p)

    @time begin
        counter = 0
        while counter < 100
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

end #module
