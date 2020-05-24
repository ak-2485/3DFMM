"""
A few ensembles for testing.
"""
module EnsembleTests

using Random, Distributions

export particledist2, particledist1, particledist3, rndmparticledist3
export particledist1ex, clusterdist

function clusterdist()
    """
    """
    particles = Dict([(1,((-4.5,-4.5,-4.5),3.0)),(2,((4.5,-4.5,-4.5),4.0)),
    (3,((-4.5,4.5,-4.5),3.0)),(4,((-4.5,-4.5,4.5),2.0)), (5,((4.5,4.5,4.5),2.0)),
    (6,((4.5,4.5,-4.5),3.0))])

    minbound = [-10.0,-10.0,-10.0]
    maxbound = [10.0,10.0,10.0]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}()
    for (key,value) in particles
        coords = value[1]
        q = value[2]
        particledict[key] = (coords,q,0.0)
    end

    return particles, particledict, minbound, maxbound
end

function particledist1()
    """
    """
    a = 2.0
    b = 2.0
    c = 2.0
    coordsx = [x for x in -19.5:a:19.5 for y in -19.5:b:19.5 for z in -19.5:c:19.5]
    coordsy = [y for x in -19.5:a:19.5 for y in -19.5:b:19.5 for z in -19.5:c:19.5]
    coordsz = [z for x in -19.5:a:19.5 for y in -19.5:b:19.5 for z in -19.5:c:19.5]

    d = 20.0
    minbound = [-d,-d,-d]
    maxbound = [d,d,d]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()

    sz=length(coordsx)
    coords = zeros(sz,3)
    coords[:,1] = coordsx
    coords[:,2] = coordsy
    coords[:,3] = coordsz

    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = 0.25*(-1)^i *(x+y)*z
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, particledict, minbound, maxbound
end

function rndmparticledist3(num::Int64)
    """
    """
    d = 10.0
    minbound = [-d,-d,-d]
    maxbound = [d,d,d]

    coords = rand(Uniform(-d,d), num,3)
    sz=size(coords)[1]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()
    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = rand(-2:0.2:2)
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, particledict, minbound, maxbound
end

end #module
