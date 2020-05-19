
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

minbound = [0.0,0.0,0.0]
maxbound = [1.0,1.0,1.0]

function mcoefTest()
    """
    Test base case multipole expansions: for box 0, entire domain.
    """
    println("Testing multipole expansion")
    p = 2
    nmax = 8 # just has to be greater than particle list tested to get single box
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    box = grid.boxes[0]

    Mnm = zeros(ComplexF64,(p+1)^2) #initialize; contains zeros at first

    # Calculate multipoles directly
    for n = 0:p
        for m = -n:n
            for (partid, vals) in particles
                q = vals[2]
                ρ, θ, ϕ = carttosphere(partid,box,grid)
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
    """
    Test the multipole expansion translation operator by checking that
    the evaluation of the two expansions are equal at a relevant point.
    """
    p = 2
    nmax = 2
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    # set the multipole expansion coefficients for  a leaf box
    child = grid.boxes[2]
    @assert !BoxStruct.haschildren(child)
    GridStruct.mcoef!(child,p,grid)
    parent = grid.boxes[child.parent]
    Onm = child.multipole_coef

    function expected(parent,child,p)
        """
        Calculate the expected translation coefficients
        """
        ρ,θ,ϕ = BoxStruct.spherecenter(child,parent.center)
        Ynm = spherical_harmonics(θ,ϕ)
        Mjk = zeros(ComplexF64,length(Onm))
        for j = 0:p
            for k = -j:j
                ind = lm_to_spherical_harmonic_index(j,k)
                for n = 0:j
                    for m = max(k+n-j,-n):min(k+j-n,n)
                        ind2 = lm_to_spherical_harmonic_index(j-n,k-m)
                        ind3 = lm_to_spherical_harmonic_index(n,-m)
                        num = Onm[ind2] * (-1)^m * BoxStruct.anms(n,m) * BoxStruct.anms(j-n,k-m) * ρ^n * Ynm[ind3]
                        Mjk[ind] += num/BoxStruct.anms(j,k)
                    end
                end
            end
        end
        return Mjk
    end
    e1 = expected(parent,child,p)

    function returned(parent,child,grid,p)
        """
        Calculate the translation coefficients from function
        """
        Mjk = BoxStruct.mcoeftrans(child,parent,p)
        return Mjk
    end
    e2 = returned(parent,child,grid,p)

    println("Testing multipole translation coefficients")

    num = 0
    den = 0
    for i in length(e1)
        num += abs(e1[i] - e2[i])^2
        den += abs(e1[i])^2
    end
    E = sqrt(num/den)
    @test E < 10^-6

    println("Testing multipole translation evaluations")

    parent.multipole_coef = e1
    point = parent.center .+ (parent.max_bound-parent.min_bound)

    Mchild = child.multipole_coef
    ψ1 = 0
    ρ1, θ1, ϕ1 = carttospherepoint(point, child.center, grid)
    Ynm = spherical_harmonics(p, θ1, ϕ1)
    for n = 0:p
        for m = -n:n
            ind1 = lm_to_spherical_harmonic_index(m,n)
            ψ1 += Mchild[ind1] * Ynm[ind1] * ρ1^-(n+1)
        end
    end

    Mparent = parent.multipole_coef
    ψ2 = 0
    ρ2, θ2, ϕ2 = carttospherepoint(point, parent.center, grid)
    Ynm = spherical_harmonics(p, θ2, ϕ2)
    for n = 0:p
        for m = -n:n
            ind1 = lm_to_spherical_harmonic_index(m,n)
            ψ2 += Mparent[ind1] * Ynm[ind1] * ρ2^-(n+1)
        end
    end

    println(" child expansion eval: ", ψ1)
    println(" parent expansion eval: ", ψ2)

    num = abs(ψ1-ψ2)^2
    den = abs(ψ1)^2
    Er = sqrt(num/den)
    @test Er < 10^-6


end

mtranslationTest()

function mtolconversionTest()
    p = 2
    nmax = 3
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
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
