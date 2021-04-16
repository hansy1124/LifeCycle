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
