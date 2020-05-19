module 3DFMM

function 3dfmm(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
    p::Int64, nmin::Int64)
     """
     """
     particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}()
     for (key,value) in particles
         coords = value[1]
         q = value[2]
         particledict[key] = (coords,q,0.0)
     end
     grid = step0(particledict,nmin)
     step1(grid,p)
     step2(grid,p)
     step3(grid,p)
     step4(grid,p)
     step5(grid,p)
     step6(grid,p)
     step7(grid,p)
     step8(grid,p)
end

end
