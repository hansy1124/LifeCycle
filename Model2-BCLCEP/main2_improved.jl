#Packages
using Plots
using Interpolations
using BenchmarkTools
using ProgressMeter
using Parameters
using StaticArrays

#Performance Packages
using Profile
using Traceur

#You can only control MyParam1 and MyParam2 for other experiments
#Then, accordingly, all codes will change

#1. Parameter Section
#1-1. Parameter for Economic Problem
MyParam1 = @with_kw (T = 100,
             β = 0.98,
             minc = 0.1,
             r = 1.02,
             σ = 2.0)
#-Can construct param while chaning only some variables
MyParam1c = MyParam1()
typeof(MyParam1c)
#1-2. Parameter for Numerical Computation
MyParam2 = @with_kw (minx=0.,
                     maxx=1000000.,
                     nX=2000,
                     nA=1000)
#-Can construct param while chaning only some variables
MyParam2c = MyParam2()
typeof(MyParam2c)

#2. (State) X-grid Generator
function gridgenerator(; minx=MyParam2c.minx,
                         maxx=MyParam2c.maxx,
                         nX = MyParam2c.nX)
    return range(minx, stop=maxx, length=nX)
end
#-Can construct gridX while changing only some variables
const gridX=gridgenerator()

#3. Life-cycle Earning Generator
lpf = [-(x-45.)^2 + 5000. for x in 1:MyParam2c.T ]
lpf[60:end] .= lpf[59]*0.6
append!(lpf,0)
typeof(lpf)
lpf = SVector{MyParam1c.T+1}(lpf) # convert it into staticarrays
#plot(lpf, title="Life-time Labor Income", label="Labor Income")

#4. current utility function
function uc(c::Vector{Float64}; σ=2.0)
    return (c.^(1-σ))./(1-σ)
end
@code_warntype uc.(maxa .- gridA)
@code_warntype uc(maxa .- gridA)
typeof(maxa .- gridA)

############################################################################
#cf) Important - where I found a critical error for this code
# - @code_warntype is really good tool to use!

#1) Usage of uc()
#Let's assume one hypothetical situation
age=9
maxa=gridX[10]
gridA=zeros(MyParam2c.nA)
gridA .= range(age == MyParam1c.T ? 0. : -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r), stop=maxa, length=MyParam2c.nA)

plot(uc.(maxa.-gridA))
plot!(uc(maxa.-gridA))
uc.(maxa.-gridA) == uc(maxa.-gridA)
# Both of the operations above have exactly same result.
# However, the types are different
@code_warntype uc.(maxa .- gridA) #return "AbstractVector"
@code_warntype uc(maxa .- gridA) #return "Vector{Float64}"
# This slowers the speed of my code significantly..!

#2) Usage of state_grid_iter! or argmax
V = zeros(MyParam1c.T+1,MyParam2c.nX)
A = zeros(MyParam1c.T,MyParam2c.nX)
age = 99
IV1 = LinearInterpolation(gridX, V[age+1, :], extrapolation_bc=Line())
@code_warntype state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, lpf, gridX) # return Tuple{Matrix{Float64}, Matrix{Float64}} - good!

@code_warntype argmax(uc(maxa .- gridA)) #return Int64 - good!

@code_warntype solver(MyParam1c, MyParam2c, lpf, gridX)
############################################################################

#5. state loop - policy loop by broadcasting
function state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, lpf, gridX)
    # For each x state, generate the choice var grid - (borrowing constraint)
    gridA=zeros(MyParam2c.nA)
    value1=zeros(MyParam2c.nA)
    for (i, x) ∈ enumerate(gridX)
        maxa = x
#        if age == MyParam1c.T
#            mina = 0.
#        else
#            mina = -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r)
#            #mina = 0
#        end
        gridA .= range(age == MyParam1c.T ? 0. : -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r), stop=maxa, length=MyParam2c.nA)
        value1 .= uc(x .- gridA) + MyParam1c.β*IV1(gridA.*MyParam1c.r.+lpf[age+1])
        p = argmax(value1)
        V[age,i], A[age,i] = value1[p], gridA[p]
    end
    return V, A
end

#6. age loop - final solver
function solver(MyParam1c, MyParam2c, lpf, gridX)
    V = zeros(MyParam1c.T+1,MyParam2c.nX)
    A = zeros(MyParam1c.T,MyParam2c.nX)
    for age in MyParam1c.T:-1:1 #@showprogress
        IV1 = LinearInterpolation(gridX, V[age+1, :], extrapolation_bc=Line())
        #IV1 = CubicSplineInterpolation(gridX, V[age+1, :], extrapolation_bc=Throw())
        V, A = state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, lpf, gridX)
    end
    return V, A
end

#V22, A22 = solver()
@btime V22, A22 = solver(MyParam1c, MyParam2c, lpf, gridX)

#@trace solver(MyParam1c, MyParam2c, lpf, gridX)
#@profile solver(MyParam1c, MyParam2c, lpf, gridX)
#Profile.print()
############################################################################

#7. Generate the Result

#7-1. Value as a function of state var at every age.
#cf) It is werid but 24th age V value doesn't make sense so
V22[24,1] = V22[23,1] # As a patch-work
plot()
for i ∈ 1:1:MyParam1c.T
    plot!(gridX, V22[i,:], label="$i Value", title="Value at Each Age")
end
plot!()

#7-2. Saving as a function of state var at every age.
plot()
for i ∈ 1:1:MyParam1c.T
    plot!(gridX, A22[i,:], label="$i Saving", title="Saving at Each Age")
end
plot!()

#7-3. Saving rate as a function of state var at every age [Negative means borrowing]
plot()
for i ∈ 1:1:MyParam1c.T
    plot!(gridX, A22[i,:]./gridX, label="$i Saving", title="Saving Rate at Each Age")
end
plot!()

#8. Life Simulation
xinit = 1000.0 #initial x
value_life = zeros(MyParam1c.T)
wealth_life = zeros(MyParam1c.T)
cons_life = zeros(MyParam1c.T)
savs_life = zeros(MyParam1c.T)

#8-1. Generate the value, consumption, saving vector
xstate = xinit
for i ∈ 1:1:MyParam1c.T
    ValueInterp=LinearInterpolation(gridX, V22[i, :], extrapolation_bc=Line())
    SavingInterp=LinearInterpolation(gridX, A22[i, :], extrapolation_bc=Line())
    wealth_life[i] = xstate
    value_life[i] = ValueInterp(xstate)
    savs_life[i] = SavingInterp(xstate)
    cons_life[i] = xstate - savs_life[i]
    xstate = savs_life[i]*MyParam1c.r + lpf[i+1]
end

#8-2. Plot the life-cycle pattern
fig1=plot(wealth_life[:], label="Life-time Wealth", title="Initial x is $xinit", fmt=:png)
plot!(savs_life, label="Saving")
plot!(cons_life, label="Consumption")
plot!(lpf, label="Labor Income")

#8-3. Save the plot
#cd(raw"C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP")
savefig(fig1,"x=$xinit.png")
