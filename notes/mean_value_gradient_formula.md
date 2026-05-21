# Mean-Value Directed Gradient in the Poincare Disk

This note documents the formula used by `directed_edge_weights.py` and
`poincare_harmonic_map.py` for directed mean-value edge weights.

## Poincare Disk Model

All vertices lie in the unit disk

$$\mathbb D = \{z \in \mathbb C : |z| < 1\}$$

The conformal factor is

$$\lambda_z = \frac{2}{1-|z|^2}$$

The Mobius addition used by the implementation is

$$x \oplus y = \frac{x+y}{1+\overline{x}y}$$

The hyperbolic distance is

$$d_{\mathbb H}(x,y) = 2\tanh^{-1}\left|(-x)\oplus y\right|$$

## Log Map

Let

$$\delta = (-x)\oplus y = \frac{y-x}{1-\overline{x}y}$$

If $|\delta|$ is nonzero, the implementation uses

$$\log_x(y) = \frac{2}{\lambda_x}\frac{\tanh^{-1}|\delta|}{|\delta|}\delta$$

Since $\lambda_x = 2/(1-|x|^2)$, this is equivalent to

$$\log_x(y) = (1-|x|^2)\frac{\tanh^{-1}|\delta|}{|\delta|}\delta$$

This tangent vector has hyperbolic norm

$$\|\log_x(y)\|_{\mathbb H,x} = d_{\mathbb H}(x,y)$$

## Exp Map

For a tangent vector $v \in T_x\mathbb D$, the implementation uses

$$\exp_x(v) = x \oplus \left(\frac{\tanh(\lambda_x |v|/2)}{|v|}v\right)$$

with $\exp_x(0)=x$.

Equivalently,

$$\exp_x(v) = \frac{x + \eta}{1+\overline{x}\eta}$$

where

$$\eta = \frac{\tanh(\lambda_x |v|/2)}{|v|}v$$

## Directed Mean-Value Weights

For a vertex $p$ and cyclic adjacent vertices $q_1,\dots,q_n$, move $p$ to the
origin:

$$u_i = \phi_p(q_i) = \frac{q_i-p}{1-\overline{p}q_i}$$

Let

$$r_i = d_{\mathbb H}(p,q_i) = 2\tanh^{-1}|u_i|$$

Let $\alpha_i$ be the angle between the geodesic directions from $p$ to $q_i$
and from $p$ to $q_{i+1}$. The implementation uses the smaller unoriented angle
between these two directions, folded into $[0,\pi]$.

The unnormalized directed mean-value weight centered at $p$ for neighbor $q_i$
is

$$w_i = \frac{\tan(\alpha_{i-1}/2)+\tan(\alpha_i/2)}{\sinh r_i}$$

For normalized weights, the implementation uses

$$\lambda_i = \frac{w_i}{\sum_j w_j}$$

The same gradient formula below applies to both $w_i$ and $\lambda_i$; only the
scale of the weights changes.

## Mean-Value Gradient Contribution

For an edge $(i,j)$, the solver stores two directed weights:

$$w_{ij}$$

centered at vertex $i$, and

$$w_{ji}$$

centered at vertex $j$.

Let

$$d_{ij} = d_{\mathbb H}(z_i,z_j)$$

The mean-value residual uses the vector

$$\sinh(d_{ij})\,e_{ij}$$

where $e_{ij}$ is the unit tangent direction at $z_i$ pointing toward $z_j$.

But the implemented log map satisfies

$$\log_{z_i}(z_j) = d_{ij}\,e_{ij}$$

Therefore the vector used in the gradient contribution is

$$\frac{\sinh(d_{ij})}{d_{ij}}\log_{z_i}(z_j)$$

The raw directed mean-value gradient update is

$$\operatorname{grad}_i \mathrel{-}= w_{ij}\frac{\sinh(d_{ij})}{d_{ij}}\log_{z_i}(z_j)$$

and

$$\operatorname{grad}_j \mathrel{-}= w_{ji}\frac{\sinh(d_{ij})}{d_{ij}}\log_{z_j}(z_i)$$

When $d_{ij}$ is numerically zero, the implementation uses the limiting scale

$$\frac{\sinh(d_{ij})}{d_{ij}} \to 1$$

## Why the Scale Is Needed

The mean-value formula divides by $\sinh r_i$, so it cancels vectors of the form

$$\sinh(r_i)e_i$$

The solver's `log_map` returns a vector of hyperbolic length $r_i$:

$$\log_p(q_i) = r_i e_i$$

Thus, to use `log_map` while preserving the mean-value identity, the solver must
multiply by

$$\frac{\sinh r_i}{r_i}$$

This factor is independent of whether the weights are normalized.
