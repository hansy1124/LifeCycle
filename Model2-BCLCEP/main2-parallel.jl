#Packages
using Distributed # Parallel computing
addprocs(4) # I have 8 processors
#Proceed after all processors are activated - kernel may break down
print(workers()) # Number of parallel workers

@everywhere using Plots
@everywhere using Interpolations
@everywhere using BenchmarkTools
@everywhere using ProgressMeter
@everywhere using Parameters
@everywhere using SharedArrays # Parallel computing

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
@everywhere MyParam2 = @with_kw (minx=0.,
                     maxx=1000000.,
                     nX=2000,
                     nA=1000)
#-Can construct param while chaning only some variables
@everywhere MyParam2c = MyParam2()
typeof(MyParam2c)

#2. (State) X-grid Generator
@everywhere function gridgenerator(; minx=MyParam2c.minx,
                         maxx=MyParam2c.maxx,
                         nX = MyParam2c.nX)
    return collect(range(minx, stop=maxx, length=nX))
end
#-Can construct gridX while changing only some variables
@everywhere gridX=gridgenerator()

#3. Life-cycle Earning Generator
@everywhere lpf = [-(x-45.)^2 + 5000. for x in 1:100 ]
@everywhere lpf[60:end] .= lpf[59]*0.6
@everywhere append!(lpf,0)
#plot(lpf, title="Life-time Labor Income", label="Labor Income")

#4. current utility function
@everywhere function uc(c; σ=MyParam1c.σ)
    return (c.^(1-σ))./(1-σ)
end

############################################################################

#5. state loop - policy loop by broadcasting
@everywhere function state_grid_iter!(V, A, IV1, age)
    # For each x state, generate the choice var grid - (borrowing constraint)
    @sync @distributed for i ∈ eachindex(gridX)
        x = gridX[i]
        maxa = x
        if age == MyParam1c.T
            mina = 0
        else
            mina = -(lpf[age+1] - MyParam1c.minc)/(1+MyParam1c.r)
            #mina = 0
        end
        nA = 1000
        gridA = collect(range(mina, stop=maxa, length=MyParam2c.nA))
        value1 = uc.(x .- gridA) + MyParam1c.β*IV1(gridA.*MyParam1c.r.+lpf[age+1])
        p = argmax(value1)
        V[age,i], A[age,i] = value1[p], gridA[p]
    end
    return V, A
end

#6. age loop - final solver
@everywhere function solver()
    V = SharedArray{Float64}(zeros(MyParam1c.T+1,MyParam2c.nX))
    A = SharedArray{Float64}(zeros(MyParam1c.T,MyParam2c.nX))
    for age in MyParam1c.T:-1:1 #@showprogress
        IV1 = LinearInterpolation(gridX, V[age+1, :], extrapolation_bc=Line())
        #IV1 = CubicSplineInterpolation(gridX, V[age+1, :], extrapolation_bc=Throw())
        V, A = state_grid_iter!(V, A, IV1, age)
    end
    return V, A
end

#V22, A22 = solver()
@btime V22, A22 = solver()

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
cd(raw"C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP")
savefig(fig1,"x=$xinit.png")
