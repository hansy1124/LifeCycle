#Packages
using Plots
using Interpolations
using BenchmarkTools

#Parameter Section
T=100
β=0.98
minc=0.1
r=1.02

#current utility function
function uc(c; σ=2.0)
    a = (c.^(1-σ))./(1-σ)
    return a
end

#X-grid Generator
maxx = 50000.
minx = 0.
nX = 1000
gridX = range(minx, stop=maxx, length=nX)
gridX = collect(gridX)

#Life-cycle Earning Generator
lpf = [-(x-45.)^2 + 5000. for x in 1:100 ]
lpf[65:end] .= lpf[64]*0.6
append!(lpf,0)

#state loop - policy loop by broadcasting
function state_grid_iter(V, A, IV1, age)
    for (i, x) ∈ enumerate(gridX)
        maxa = x
        if age == 100
            mina = 0
        else
            mina = -(lpf[age+1] - minc)/(1+r)
            #mina = 0
        end
        nA = 1000
        gridA = range(mina, stop=maxa, length=nA)
        gridA = collect(gridA)
        value1 = uc.(x .- gridA) + β*IV1(gridA.*r.+lpf[age+1])
        p = argmax(value1)
        V[age,i], A[age,i] = value1[p], gridA[p]
    end
    return V, A
end

function solver()
    V = zeros(T+1,nX)
    A = zeros(T,nX)
    for age in T:-1:1
        IV1 = LinearInterpolation(gridX, V[age+1, :], extrapolation_bc=Line())
        V, A = state_grid_iter(V, A, IV1, age)
    end
    return V, A
end

V22, A22 = solver()


plot(lpf)

plot()
for i ∈ 1:1:T
    plot!(gridX, V22[i,:], label="$i Value")
end
plot!()

plot()
for i ∈ 1:1:T
    plot!(gridX, A22[i,:], label="$i Saving")
end
plot!()

#Life - Simulation
xstate = 1000.0
value_life = zeros(T)
wealth_life = zeros(T)
cons_life = zeros(T)
savs_life = zeros(T)
for i ∈ 1:1:T
    ValueInterp=LinearInterpolation(gridX, V22[i, :], extrapolation_bc=Line())
    SavingInterp=LinearInterpolation(gridX, A22[i, :], extrapolation_bc=Line())
    wealth_life[i] = xstate
    value_life[i] = ValueInterp(xstate)
    savs_life[i] = SavingInterp(xstate)
    cons_life[i] = xstate - savs_life[i]
    xstate = savs_life[i]*r + lpf[i+1]
end

fig1=plot(wealth_life[:], label="Life-time Wealth", title="Initial x is 1000", fmt=:png)
plot!(savs_life, label="Saving")
plot!(cons_life, label="Consumption")
plot!(lpf, label="Labor Income")

#cd(raw"C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP")
savefig(fig1,"x=1000.png")
