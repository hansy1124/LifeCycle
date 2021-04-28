#Note) GR Plots QuantEcon packages are conflicting to each other
#Build GR first (and do not use it), then declare "using Plots -> using QuantEcon"

#Library
using Distributed # Parallel computing
addprocs(8) # I have 8 processors
#Proceed after all processors are activated - kernel may break down
print(workers()) # Number of parallel workers

#@everywhere using GR
@everywhere using Plots
@everywhere using Interpolations
@everywhere using BenchmarkTools
@everywhere using ProgressMeter
@everywhere using Parameters
@everywhere using QuantEcon
@everywhere using FastGaussQuadrature
@everywhere using LinearAlgebra
@everywhere using SharedArrays

include("parameter-parallel.jl")
include("model-parallel.jl")
#Parameters fro economic problem
@everywhere MyParam1c = MyParam1()

#Parameters for numerical sol
@everywhere MyParam2c = MyParam2(nX=500, nY = 100, nA = 500)

#Labor Income Shocks
@everywhere v, w = normal_gausshermite(MyParam1c.μ_v, MyParam1c.σ_v, MyParam1c.nv)

#Labor Income Growth Rate
#cf) gbar[age] : growth rate from Y_age to Y_age+1
@everywhere gbar, lpf = generate_gbar(MyParam1c.T, MyParam1c.reT)
p1 = plot(gbar, title="Labor Income Growth Rate", label="", xlabel="Age")
p2 = plot(lpf, title="Labor Income Profile", label="", xlabel="Age")
plot(p1, p2, layout=2)

#Generating Grids
#cf) Generate the grid for X
@everywhere gridX = range(MyParam2c.minx,
              stop=MyParam2c.maxx,
              length=MyParam2c.nX)

#cf) Generate the max min grid point for Y - Age specific grid will be generated
@everywhere VmaxY, VminY = Ygrid_life(1000., lpf, v, MyParam1c)
plot(lpf, title="Life Cycle Earning Profile", label="Expected")
plot!(VmaxY, label="Max-Y Path")
plot!(VminY, label="Min-Y Path")

#Given an age
@everywhere function state_grid_iter!(V, A, IV1, age,
                          MyParam1c, MyParam2c, gbar,
                          gridX, gridY, gridT)
    #Grid for optimal policy
    gridA = zeros(MyParam2c.nA)
    #Container for values
    value1 = zeros(MyParam2c.nA)

    #Loop over states
    @sync @distributed for i in 1:1:size(gridT)[1]
        xind = Int64(gridT[i,:][1]) # x grid point
        yind = Int64(gridT[i,:][2]) # y grid point
        x = gridX[xind]
        y = gridY[yind]

        #Generate a choice grid for age/state specific
        maxa = x
        #cf) gbar[age] : growth rate from y_age to y_{age+1}
        gridA .= collect(range(age == MyParam1c.T ? 0. : -(y*exp(gbar[age]+v[1]) - MyParam1c.minc)/MyParam1c.r,
        stop=maxa,
        length=MyParam2c.nA))
        #Even in the worst case, household must be able to pay back the debt and cosnume the minimum level of consumption
        #gridA .= range(0., stop=maxa, length=MyParam2c.nA) #for no borrowing

        #Given a set of state variables / grid for choice
        #- loop across the choice var
        for j in eachindex(gridA)
            value1[j] = uc(x - gridA[j]) + MyParam1c.β*dot(
            IV1.(MyParam1c.r*gridA[j].+y.*exp.(gbar[age].+v), #1st arg
            y.*exp.(gbar[age].+v)), #2nd arg
            w) #prob dist
            #cf) for IV1. (Interpolation function)
            #1st argument - MyParam1c.r*gridA[i].+y.*exp.(gbar[age].+v)
            #2nd argument - y.*exp.(gbar[age].+v)
            #cf) for dot() - Expectation
            #1st argument IV1.()
            #2nd argument w
        end
        p = argmax(value1)
        V[age,xind,yind], A[age,xind,yind] = value1[p], gridA[p]
    end
    return V, A
end

#6. Age loop - final solver

@everywhere function solver(MyParam1c, MyParam2c, gbar, gridX, VmaxY, VminY)
    V = SharedArray{Float64}(zeros(MyParam1c.T+1, MyParam2c.nX, MyParam2c.nY))
    A = SharedArray{Float64}(zeros(MyParam1c.T, MyParam2c.nX, MyParam2c.nY))
    gridT=gridmake(range(1, stop=MyParam2c.nX, length=MyParam2c.nX),
                   range(1, stop=MyParam2c.nY, length=MyParam2c.nY))
    #Grid container for each age
    #gridTC = [Matrix{Float64}(undef,dim1,dim2) for i in 1:MyParam1c.T]
    @showprogress for age in MyParam1c.T:-1:1
        #Make a age-specific one-period future Y-grid for interpolation
        gridY1=collect(range(VminY[age+1],
                       stop=VmaxY[age+1]+0.1, # for the last period
                       length=MyParam2c.nY))
        #Make a age-specific Y-grid
        gridY=collect(range(VminY[age],
                      stop=VmaxY[age],
                      length=MyParam2c.nY))
        #Interpolate
        V2 = [V[age+1,x,y] for x in eachindex(gridX), y in eachindex(gridY1)]
        IV1 = LinearInterpolation((gridX, gridY1),
                                   V2,
                                   extrapolation_bc=Line())
        #Store a age specific total state space
        V, A = state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, gbar, gridX, gridY, gridT)
    end
    return V, A
end

V, A = solver(MyParam1c, MyParam2c, gbar, gridX, VmaxY, VminY)

##############################################
############### 6. Show Result ###############
##############################################

#1) Should be all zero (Before you die, You consume all)
print(A[MyParam1c.T,2:40,:])

#2) Checking for each age
ageC = MyParam1c.T-2
#1-Value for each state space across X and Y
print(V[ageC,:,50])

#1-1. 3D
#gridX
y1 = range(VminY[ageC], stop=VmaxY[ageC], length=MyParam2c.nY)
p = plot(y1, gridX, V[ageC,:,:], st=:surface, camera=(-50,50), xlabel="Income", ylabel="Wealth", zlabel="Value", title="Value at Age $ageC")
plot!(p[1], camera = (25,50))

#1-2. By Age
plot()
for i in [1,20,50,70,100]
    plot!(V[ageC,2:40,Int64(i)], label=string(i))
end
plot!()

#2-Optimal Policy for Saving across X and Y
print(A[ageC,:,50])

#2-1. 3D
#gridX
#gridX
y1 = range(VminY[ageC], stop=VmaxY[ageC], length=MyParam2c.nY)
p = plot(y1, gridX, A[ageC,:,:], st=:surface, camera=(-50,50), xlabel="Income", ylabel="Wealth", zlabel="Value", title="Saving at Age $ageC")
plot!(p[1], camera = (30,50))

#2-2. By Age
plot()
for i in [1,20,50,70,100]
    plot!(A[ageC,2:50,Int64(i)], label=string(i))
end
plot!()
#########################################################################
