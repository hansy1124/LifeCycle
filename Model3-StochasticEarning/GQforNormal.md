## Gaussian Quadrature for Normal Distribution

Following are the sources that I checked
- https://jblevins.org/notes/quadrature
- https://www.wouterdenhaan.com/numerical/integrationslides.pdf - 24th slide
- http://fmwww.bc.edu/ec-p/software/Miranda/chapt6.pdf - 5th page

#### Problem

Let's say I want to calculate the below.

$$
\begin{align*}
E[h(y)] \quad with \quad y \sim N(μ, σ^2)
\end{align*}
$$

Here, to properly calculate this, I need to do the change of variables (which I will forget so I skip)

#### Procedure

1) Obtain n Gauss-hermite quadrature weights and nodes using a numerical algorithm
2) Calculate the approximation by
$$
\begin{align*}
E[h(y)] \sim \sum_{i=1}^n {1 \over √π} w_i^{GH} h(√2 σ ζ_i^{GH}+μ)
\end{align*}
$$
3) So, the proper discrete distribution will be, for i,

Node: $√2 σ \zeta_i^{GH} + \mu$ & Prob: $w_i^{GH}$

#### My function: normal_gausshermite(μ, σ, n)
This functions give you n grid gauss_hermite quadrature for normal distribution with mean μ and std σ. 
