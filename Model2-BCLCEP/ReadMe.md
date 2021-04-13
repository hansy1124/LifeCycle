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
