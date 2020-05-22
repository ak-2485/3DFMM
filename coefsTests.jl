
using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using EnsembleTests
using BoxStruct
#using PlotGrid
import FMMSteps
import FMM3D


function mcoefTest()
    """
    Test base case multipole expansions: for box 0, entire domain.
    """
    println("Testing multipole expansion")
    p = 2
    nmax = 8 # just has to be greater than particle list tested to get single box
    particles, minbound, maxbound = particledist1()
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

#mcoefTest()

function mtranslationTest1()
    """
    Test the multipole expansion translation operator by checking that
    the evaluation of the two expansions are equal at a relevant point.
    """
    p = 1
    nmax = 2
    particles, minbound, maxbound = particledist1()
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    #plotgrid(grid)
    child = grid.boxes[2]
    @assert child.numparticles==1
    a = copy(child.particles)
    s = pop!(a)
    println(grid.particles[s])
    println(child.center)
    @assert !BoxStruct.haschildren(child)
    # set the multipole expansion coefficients for  a leaf box
    GridStruct.mcoef!(child,p,grid)
    parent = grid.boxes[child.parent]
    println(parent.center)
    Onm = child.multipole_coef

    function expected(parent,child,p)
        """
        Calculate the expected translation coefficients
        """
        ρ,θ,ϕ = BoxStruct.spherecenter(child,parent.center)
        Ynm = spherical_harmonics(p,θ,ϕ)
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

    parent.multipole_coef = e2
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
    println("Error is: ", Er)
end

#mtranslationTest1()

function mtranslationTest2()
    """
    Compute the multipole coefficients of the grid as one box with no refinement
    and then again as a pass from 8 children to parent. Compare at relevant point.
    """
    p = 1
    nmax = 6
    particles, minbound, maxbound = particledist1()
    grid1 = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    FMM3D.step1(grid1,p)
    println("new grid")
    p = 1
    nmax = 2
    grid2 = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    FMM3D.step1(grid2,p)
    FMM3D.step2(grid2,p)

    parent1 = grid1.boxes[0]
    parent2 = grid2.boxes[0]

    point = parent1.center .+ (parent1.max_bound-parent1.min_bound)
    ψ1 = 0
    ψ2 = 0
    ρ1, θ1, ϕ1 = carttospherepoint(point, parent1.center, grid1)
    ρ2, θ2, ϕ2 = carttospherepoint(point, parent2.center, grid2)
    Ynm1 = spherical_harmonics(p, θ1, ϕ1)
    Ynm2 = spherical_harmonics(p, θ2, ϕ2)
    Mparent1 = parent1.multipole_coef
    Mparent2 = parent2.multipole_coef
    for n = 0:p
        for m = -n:n
            ind1 = lm_to_spherical_harmonic_index(m,n)
            ψ1 += Mparent1[ind1] * Ynm1[ind1] * ρ1^-(n+1)
            ψ2 += Mparent2[ind1] * Ynm2[ind1] * ρ2^-(n+1)
        end
    end

    num = abs(ψ1-ψ2)^2
    den = abs(ψ1)^2
    Er = sqrt(num/den)
    println("Error is: ", Er)

end

#mtranslationTest2()

function mtranslationTest3()

    particles, particledict, minbound, maxbound = particledist3()

    p = 1
    nmax = 50
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    #plotgrid(grid)
    FMM3D.step1(grid,p)
    parent = grid.boxes[0]

    point = parent.center .+ (parent.max_bound-parent.min_bound)

    for childid in parent.children
        child = grid.boxes[childid]
        Mjk = BoxStruct.mcoeftrans(child,parent,p)
        Mchild = child.multipole_coef
        ψ1 = 0
        ρ1, θ1, ϕ1 = carttospherepoint(point, parent.center, grid)
        Ynm = spherical_harmonics(p, θ1, ϕ1)
        for n = 0:p
            for m = -n:n
                ind1 = lm_to_spherical_harmonic_index(m,n)
                ψ1 += Mchild[ind1] * Ynm[ind1] * ρ1^-(n+1)
            end
        end

        Mparent = Mjk
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

        num = abs(ψ1-ψ2)
        den = abs(ψ1)
        Er = sqrt(num/den)
        println("Error is: ", Er)
    end


end #test
#mtranslationTest3()

function mtomTest4()
    """
    Test the M2M operator by first comparing the moments of a parent
    box to the moments of the same parent box after translation from
    children.
    """
    println("Testing m2m: parent <--- children")
    particles, minbound, maxbound = clusterdist()
    p = 2
    println("p = ", p)
    nmax = length(particles)
    grid1 = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    # set the moments of the parent box when it is childless
    FMM3D.step1(grid1,p)
    # on a new grid, set the moments of the children then pass them to parent
    nmax = 2
    grid2 = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    FMM3D.step1(grid2,p)
    FMM3D.step2(grid2,p)
    # compare moments of the parent
    parent1 = grid1.boxes[0]
    parent2 = grid2.boxes[0]
    Mparent1 = parent1.multipole_coef
    Mparent2 = parent2.multipole_coef
    println(Mparent1)
    println(Mparent2)
    num = 0
    den = 0
    for i = 1:length(Mparent1)
        Mnmexpected = Mparent1[i]
        Mnm = Mparent2[i]
        num += abs(Mnmexpected - Mnm)^2
        den += abs(Mnmexpected)^2
    end
    E = sqrt(num/den)
    println("Error is: ", E)

end

mtomTest4()

function ltolTest1()
    """
    Test the local to L2L operator by comparing the lcoefs of the parent due to
    a colleague translated to a child and the lcoefs of the child due to the
    same colleague of the parent. These values are equal for p=0.
    """
    println("Testing l2l: parent --> child due to parent colleague")
    particles,_,minbound, maxbound = particledist1ex()
    # This test requires p = 0 !
    p = 0
    nmax = 20
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    # pick a parent box from the grid
    parent = grid.boxes[1]
    # get a colleague
    colparent = grid.boxes[19]
    @assert parent.level == colparent.level
    # set the local coefs of the parent due to colleague's particles
    lcoef!(parent,colparent,p,grid)
    # get array of coefs that would pass to the child
    child = grid.boxes[pop!(parent.children)]
    Ljk = lcoeftrans(parent, child, p)
    # set the local coefs of the child due to colleague
    lcoef!(child,colparent,p,grid)
    println(Ljk)
    println(child.local_coef)

    num = 0
    den = 0
    for i = 1:length(Ljk)
        Ljkexpected = child.local_coef[i]
        num += abs(Ljkexpected - Ljk[i])^2
        den += abs(Ljkexpected)^2
    end

    E = sqrt(num/den)
    println("Error is: ", E)

end

#ltolTest1()
