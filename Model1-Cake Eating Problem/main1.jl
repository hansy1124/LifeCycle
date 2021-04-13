"""
Cake Eating Problem

- Intertemporal problem of consumer living for T periods and endowed with initial wealth a_1 in period t=1
- To allocate the consumption of this wealth over her T periods of life in order to maximise her lifetime wellbeing
- Consumption is divisible: a continuous decision variable
- Any remaining wealth in period t is productive, generating k(a) units of wealth to consume in the future
- No outstanding debts are allowed at the end of life
- And ay remaining wealth at the end of life is of no value
"""

# type Ctrl+D to shut down the kernel and restart

using Plots

#########################################################################

# From Euler Equation, we know that the analytical solution is
function cakeeating(T::Int64, a1::Float64, β::Float64, R::Float64, γ::Float64)
    α = β^(1/γ)*R^((1-γ)/γ)
    c = zeros(T)
    for i in eachindex(c)
        c[i] = (1-α)*a1*(β*R)^((i-1)/γ)/(1-α^T)
    end
    return c
end

#########################################################################

# Solve one specific case
T=100
a1=100.0
β=0.98
R=1.04
γ=1.5
sol1 = cakeeating(T, a1, β, R, γ)
plot(1:T, sol1, label="C* w/ γ=$γ", title="Cake Eating Prob β=$β, and R=$R", xlabel="Age", ylabel="Consumption Amount")
plot()

#########################################################################

# Solve for several R valeus given β = 0.98
R1 = [i for i in 1.01:0.01:1.03]
for j in eachindex(R1)
    sol1 = cakeeating(T, a1, β, R1[j], γ)
    R2=R1[j]
    plot!(1:T, sol1, label="C* w/ R=$R2", title="Cake Eating Problem given β=$β and γ=$γ", xlabel="Age", ylabel="Consumption Amount")
end
plot!()
plot()

#########################################################################

# Solve for several γ (CRR Coefficient and EIS) values given β=0.98 and R=?
β=0.98
R=1.03
γ1 = [2.0i for i in 1:5]
@show γ1
for γ in γ1
    sol1 = consumption(T, a1, β, R, γ)
    plot!(1:T, sol1, label="C* w/ γ=$γ",xaxis=[0,100])
end
plot!()
plot()

#-> We can see that γ control the eagerness to smooth the consumption

#########################################################################

# In terms of optimal consumption and saving policy, since "a" is a state variable, I can plot optimal policy as a function of "a".

# By the formula, ..
t = [t for t in 1.0:99.00] # remaining age
α = β^(1/γ)*R^((1-γ)/γ)
i1 = similar(t)
i1 .= 1
@. optcshare = (i1-α)/(i1-α^t)
optashare = 1.0.-optcshare

plot(t, optcshare, label="consumption share c/a", title="C/S share as f(remaining age)", xlabel="Remaining Ages")
plot!(t, optashare, label="saving share (a-c)/a")
plot()

# As you see the further future, you consume less portion of your current asset now

#########################################################################
