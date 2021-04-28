#Packages
using Distributed # Parallel computing
addprocs(8) # I have 8 processors
#Proceed after all processors are activated - kernel may break down
print(workers()) # Number of parallel workers

@everywhere using Plots
@everywhere using Interpolations
@everywhere using BenchmarkTools
@everywhere using Parameters
@everywhere using ProgressMeter
@everywhere using StaticArrays
@everywhere using SharedArrays

#Performance Packages
using Profile
using Traceur

#You can only control MyParam1 and MyParam2 for other experiments
#Then, accordingly, all codes will change

#1. Parameter Section
#1-1. Parameter for Economic Problem
@everywhere MyParam1 = @with_kw (T = 100,
             β = 0.98,
             minc = 0.1,
             r = 1.02,
             σ = 2.0)
#-Can construct param while chaning only some variables
@everywhere MyParam1c = MyParam1()
typeof(MyParam1c)
#1-2. Parameter for Numerical Computation
@everywhere MyParam2 = @with_kw (minx=0.1,
                     maxx=50000.,
                     nX=1000,
                     nA=1000)
#-Can construct param while chaning only some variables
@everywhere MyParam2c = MyParam2()
typeof(MyParam2c)

#2. (State) X-grid Generator
@everywhere function gridgenerator(; minx=MyParam2c.minx,
                         maxx=MyParam2c.maxx,
                         nX = MyParam2c.nX)
    return range(minx, stop=maxx, length=nX)
end
#-Can construct gridX while changing only some variables
@everywhere gridX=gridgenerator()

#3. Life-cycle Earning Generator
@everywhere lpf = [-(x-45.)^2 + 5000. for x in 1:MyParam1c.T]
@everywhere lpf[60:end] .= lpf[59]*0.6
@everywhere append!(lpf,0)
typeof(lpf)
@everywhere lpf = SVector{MyParam1c.T+1}(lpf) # convert it into staticarrays
#plot(lpf, title="Life-time Labor Income", label="Labor Income")

#4. current utility function
@everywhere function uc(c::Vector{Float64}; σ=2.0)
    return (c.^(1-σ))./(1-σ)
end

#5. state loop - policy loop by broadcasting
@everywhere function state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, lpf, gridX)
    # For each x state, generate the choice var grid - (borrowing constraint)
    gridA=zeros(MyParam2c.nA)
    value1=zeros(MyParam2c.nA)
    @sync @distributed for i ∈ eachindex(gridX)
        x = gridX[i]
        maxa = x
#        if age == MyParam1c.T
#            mina = 0.
#        else
#            mina = -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r)
#            #mina = 0
#        end
        gridA .= range(age == MyParam1c.T ? 0. : -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r), stop=maxa, length=MyParam2c.nA)
        #gridA .= range(0., stop=maxa, length=MyParam2c.nA)
        value1 .= uc(x .- gridA) + MyParam1c.β*IV1(gridA.*MyParam1c.r.+lpf[age+1])
        p = argmax(value1)
        V[age,i], A[age,i] = value1[p], gridA[p]
    end
    return V, A
end

#6. age loop - final solver
@everywhere function solver(MyParam1c, MyParam2c, lpf, gridX)
    V = SharedArray{Float64}(zeros(MyParam1c.T+1,MyParam2c.nX))
    A = SharedArray{Float64}(zeros(MyParam1c.T,MyParam2c.nX))
    for age in MyParam1c.T:-1:1 #@showprogress
        IV1 = LinearInterpolation(gridX, V[age+1, :], extrapolation_bc=Line())
        #IV1 = CubicSplineInterpolation(gridX, V[age+1, :], extrapolation_bc=Throw())
        V, A = state_grid_iter!(V, A, IV1, age, MyParam1c, MyParam2c, lpf, gridX)
    end
    return V, A
end
#@time
#@btime V22, A22 = solver(MyParam1c, MyParam2c, lpf, gridX)
V22, A22 = solver(MyParam1c, MyParam2c, lpf, gridX)

#@trace solver(MyParam1c, MyParam2c, lpf, gridX)
#@profile solver(MyParam1c, MyParam2c, lpf, gridX)
#Profile.print()

############################################################################

#7. Generate the Result

# Note on error of the solution
# By construction of the function, choice variable for consumption (c=x-a) cannot be negative. However, due to numerical unstability of the julia when we calculate really small number, sometimes the code end up having negative consumption, which pushes the value to the extreme positive given the specification of utility function. For example, if the x state variable is close to 0 or literally 0, the saving choice variable unstably become 0.0000001 (from a range whose max is 0), which results in negative consumption and positive value.
#Check where this happens ..
for i ∈ 1:100
    #if V22[i,1] < -0.5
    if V22[i,1] > 0.0
        print("Problem in V from ", i, "a")
    end
end

#As a solution,
#1) you can set starting point of x as 0.1, or 1.0, or 10. - Then, you will have no such problem as c=x-a won't be negative. In addition, as the value function is not that steep, the linear interpolation will not make a big problem. (Note, since x is the beginning wealth including the labor income, at age=100, the value will be really low because this person literally has no wealth at all - it is not an error.)

#7-1. Value as a function of state var at every age.
plot()
for i ∈ 1:1:MyParam1c.T-1
    plot!(gridX, V22[i,:], label="$i Value", title="Value at Each Age")
end
plot!()

#7-2. Value as a function of state var and age (3d)
pyplot()
x=range(1, stop=MyParam1c.T-1, length=MyParam1c.T-1) #remove at 100 as it has too negative value
y=range(1, stop=length(gridX), length=length(gridX))
f(x,y) = V22[Int64(x),Int64(y)]
plot(x,y,f,st=:surface,camera=(-30,30), title="Value across X and Age")

#7-3. Saving as a function of state var at every age.
plot()
for i ∈ 1:1:MyParam1c.T
    #plot!(gridX, A22[i,:], label="$i Saving", title="Saving at Each Age")
    plot!(gridX, A22[i,:], label="", title="Saving at Each Age")
end
plot!()

#7-4. Saving as a function of state var and age (3d)
pyplot()
x=range(1, stop=MyParam1c.T-1, length=MyParam1c.T-1) #remove at 100 as it has too negative value
y=range(1, stop=length(gridX), length=length(gridX))
f(x,y) = A22[Int64(x),Int64(y)]
plot(x,y,f,st=:surface,camera=(-30,30), title="Saving across X and Age")

#7-5. Saving rate as a function of state var and age (3d) [Negative means borrowing]
A23=similar(A22)
for i ∈ 1:1:MyParam1c.T
    A23[i,:] = A22[i,:]./gridX
end
pyplot()
x=range(1, stop=MyParam1c.T-1, length=MyParam1c.T-1) #remove at 100 as it has too negative value
y=range(1, stop=length(gridX), length=length(gridX))
f(x,y) = A23[Int64(x),Int64(y)]
plot(x,y,f,st=:surface,camera=(-30,30), title="Saving rate across X and Age")

############################################################################

#8. Life Simulation

xinit = 50000.0 #initial x
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
