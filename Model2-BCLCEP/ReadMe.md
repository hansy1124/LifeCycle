## Model 2 - Borrowing Constraint and Life-Cycle Earning Profile

This file is intended to solve a life-cycle consumption and saving problem with realistic life-cycle (deterministic) labor income profile and borrowing constraint.

Following is the model. (Note that $Y_t$ is not a state variable, since we do not have stochastic labor income profile. I assume that household knows the whole future path of their labor income.)

$$\begin{align*}
V(X_t) &= max \ u(C_t) + β V(X_{t+1}) \\
s.t. \quad X_t &= C_t +A_t \\
 X_{t+1} &= A_t (1+r) + Y_{t+1} \\
 A_t &\ge -b
\end{align*}$$

Here, I won't allow the default, which means that, the borrowing constraint will be life-cycle component. At time $t$,

$$\begin{align*}
(1+r) b \le Y_{t+1} - min [C]
\end{align*}$$

cf) If the labor income becomes stochastic, it will be messier.

I will use the grid-search method for solving this model. Grid searching model becomes better as the model becomes more complicated, though the computation time becomes longer, with kinks where the usual solvers would end up giving us local maximum.

#### Note on Parallel Computing
main2.jl - file for solving the model - Basic Julia
main2-parallel.jl - file for solving the model with "Distributed Packages."
cf) I changed how I iterate over x-state to use @distributed correctly.
cf) I changed "for (i, x) ∈ enumerate(gridX)" into
"@sync @distributed for i ∈ eachindex(gridX) \ x = gridX[i]"

#### Result on Performances (through @btime)
1) Performance w/o Parallel Computing
51.941 s (908208635 allocations: 28.70 GiB)
49.259 s (908126593 allocations: 28.70 GiB)
50.021 s (908126593 allocations: 28.70 GiB)
51.126 s (908126593 allocations: 28.70 GiB)

2) Performance w/ Parallel Computing with 8 Processors
15.959 s (215836 allocations: 40.01 MiB)
14.281 s (200246 allocations: 38.81 MiB)
17.137 s (200232 allocations: 38.79 MiB)
15.699 s (200230 allocations: 38.78 MiB)

3) Performance w/ Parallel Computing with 4 Processors
20.007 s (98361 allocations: 21.51 MiB)
19.611 s (98361 allocations: 21.51 MiB)

#### Imporvement in Performance by changing codes
I documented change from main2.jl to main2_improved.jl
1) Initial performance from main2
- 48.425 s (908126593 allocations: 28.70 GiB)
2) Put parameter tuples into function as argument rather than using parameter tuples as global variables
- 46.684s (905742585 allocations: 28.63 GiB)
- 47.311 s (905742585 allocations: 28.63 GiB)
3) In state_grid_iter!, mina becomes either Int64 (0) or Float64 (calculated by next period labor income) in the original code  So I change it into 0. so that it always becomes Float64
- 47.828 s (905744585 allocations: 28.63 GiB)
4) Remove unnecessary codes such as nA=1000, V1 = V, A1 = A
- 47.337 s (905744585 allocations: 28.63 GiB)
5) Rather than using p = argmax(value1) and value1[p], use value1[argmax(value1)]
-> Code becomes worse
- 50.515 s (905939465 allocations: 28.63 GiB)
6) Change lpf as StaticArrays
- 46.648 s (905744685 allocations: 28.63 GiB)
- 47.111 s (905744685 allocations: 28.63 GiB)
7) When update gridA and value1, replace ( .=) rather than generating new arrays (=)
- 47.715 s (904350005 allocations: 28.59 GiB)
8) Put gridX as also an argument to functions
- 47.244 s (901000505 allocations: 28.50 GiB)
9) When update the values and optimal policy containers, I use (.=) replace rather than generating new one
- 48.097 s (901000505 allocations: 28.50 GiB) - no difference
10) I tried traceur. Following is the result.
┌ Warning: uses global variable Main.MyParam1c
└ @ C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP\main2_improved.jl:53
┌ Warning: dynamic dispatch to Base.getproperty(Main.MyParam1c, σ)
└ @ C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP\main2_improved.jl:53
┌ Warning: dynamic dispatch to Main.:(var"#uc#12")(Base.getproperty(Main.MyParam1c, σ), _1, _2)
└ @ C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP\main2_improved.jl:53
┌ Warning: uc returns Any
└ @ C:\Users\Owner\Dropbox\8. Coding Space\Julia\LifeCycle\Model2-BCLCEP\main2_improved.jl:52
It helped me to understand that my uc() function returns weird AnyVector type outcome. So I changed it so that uc() returns Vector{Float64} type outcome.

11. However, above one was not enough because I used uc.() in the loop. @code_warntype actually helped me to know that uc.() returns AnyVector type still. (-> Body::AbstractVector{var"#s831"} where var"#s831") So I changed the code into just uc() which returns Vector{Float64}
- 21.514 s (1200505 allocations: 9.09 GiB)

12. I used a profiling, and it shows that interpolatio and uc are the bottlenecks for this code.
- 3803  ...ace\Julia\LifeCycle\Model2-BCLCEP\main2_improved.jl:54; uc
-  6952  @Interpolations\src\extrapolation\extrapolation.jl:49; Extrapolation

13. Using ternary operator rather than if-else statement does not change the outcome much. In addition, range and linspace are same function. So no difference.

### Doing improved version in parallel computing

1) With 4 processors
- 5.332 s (143963 allocations: 26.40 MiB)

2) With 8 processors
- 2.783 s (330212 allocations: 47.93 MiB)
