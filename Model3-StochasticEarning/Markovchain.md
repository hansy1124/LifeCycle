I am going to use Tauchen and Rouwenhorst from QuantEcon.jl (tauchen(), rouwenhorst())

For correlated random variables process, I may need to do it by myself. Below is the brief description for each method.



For other codes, ...
cf) https://github.com/sglyon/Q4Models/blob/master/julia/discretize.jl
cf) https://github.com/gcamilo/tauchenHussey

#### One dimensional Markov Chain (Tauchen 1985)

##### 1) Model
Let's assume following model.
$$y_{t} = ρ y_{t-1} + ϵ_t \qquad ϵ_t \sim N(0,\sigma^2)$$

Then, we know that, assuming stationary process ($ρ < 1$)

##### 2) Moments of the process
$$\begin{align*}
E(y) &= 0 \\
SE(y) &= \sqrt{  σ^2  \over 1-\lambda^2 }
\end{align*}$$

##### 3) Discretization
###### 3-1) Grid
Then, we need to discretize the process. Ginve the number of grid for this process $(n)$ so $[\bar{y}^1, \bar{y}^2, ... , \bar{y}^{N-1}, \bar{y}^N]$, we can do following.
1. Set $\bar{y}^N = mSE(y)$ for some constant $m$
2. Set $\bar{y}^1 = -\bar{y}^N$
3. Let the remaining be equispaced over the interval $[\bar{y}^1,\bar{y}^N]$

###### 3-2) Transition Prob Matrix
For $p_{jk} = Pr[\tilde{y}_t=\bar{y}^k|\tilde{y}_{t-1}=\bar{y}^j]$, put $w = \bar{y}^k - \bar{y}^{k-1}$. Then, we can calculate the transition probability by assuming $ρy_{t-1} + ϵ_t$ as a random variable.

$$
\begin{align*}
P_{jk} &= P [\bar{y}^k - w/2 ≤ ρ\bar{y}^j + ϵ_t ≤ \bar{y}^k + w/2] \\
       &= F({\bar{y}^k-\lambda \bar{y}^j + w/2 \over σ}) - F({\bar{y}^k-\lambda \bar{y}^j - w/2 \over σ})
\end{align*}
$$

#### Multivariate Markov Chain (Tauchen 1985)
Suppose the vector process below,

$y_t = A y_{t-1} + ϵ_t, var(ϵ_t) = \Sigma_{ϵ},$

Variance matrix is a diagonal matrix - all white noises are not correlated. The dimension of each vector is $y\ (= M × 1), A\ (= M × M), ϵ\ (= M × 1)$.

Effectively, it proceeds almost same way like above. Key things are that 1) each $y_i$ in $y$ will have each $n_i$ grid points, which makes the number of total states as $N=\Pi_{i=1}^M n_i$.

Again, the transition probability will be just the **multiplication** of probabilities of all components moving from one (their) possible realization to another (their) possible realization. (Note that each state for this vector Markov chain specifies the exact realization cases for all $y_i$ for $i=1, .., M$)

Though all white noises are correlated, the resulting $y$ process can have non-diagonal covariance variance depending on the specification of $A$.

#### One dimensional Markov Chain (Tauchen Hussey 1991)
This method extends the Tauchen 1985 by introducing "quadrature method" [https://www.youtube.com/watch?v=nQZYBWB6q_k - for Gaussian Quadrature - selecting grids and weights to linear mapping can approximate the integral well] for selecting grids. Tauchen 1985 works well for simple problem (using equispaced grids). However, it cannot handle large problems efficiently including complex dynamics especially conditional heteroskedasticity.

The grid points are determined by f(y|x) (Conditional distribution) in conjunction with a numerical quadrature rule, and this method can handle large problems and more elaborate dynamics than those of a linear VAR.

check https://github.com/gcamilo/tauchenHussey

### Cross-correlated VAR (1) Discretization (Galindev Lkhagvasuren 2010)

### Kopecky Suen 2010 - For more persistent VAR(1) process
- Rouwenhort Method is well-explained.

#### Rouwenhort 1995 Method
- Use symmetric and evenly spaced grid
- Transition matrix with some special formula with $p, q$ in 703 page 
