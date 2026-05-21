Given $p\in\mathbb D$ and cyclic adjacent points

$$
q_1,\dots,q_n\in\mathbb D,\qquad q_{n+1}=q_1,
$$

first move $p$ to the origin using the Poincare disk isometry

$$
u_i
=
\phi_p(q_i)
=
\frac{q_i-p}{1-\overline p\,q_i}.
$$

Then the hyperbolic distance from $p$ to $q_i$ is

$$
r_i
=
d_{\mathbb H}(p,q_i)
=
2\tanh^{-1}|u_i|.
$$

Let $\alpha_i$ be the angle at $p$ between the geodesics from $p$ to $q_i$ and $q_{i+1}$. Since the Poincare model is conformal,

$$
\alpha_i
=
\arg(u_{i+1})-\arg(u_i)
$$

taken as the positive cyclic angle.

The **hyperbolic mean value weight** for neighbor $q_i$ is

$$
\boxed{
w_i
=
\frac{
\tan(\alpha_{i-1}/2)+\tan(\alpha_i/2)
}{
\sinh r_i
}
}
$$

with indices taken cyclically.


The normalized mean value coordinate is

$$
\boxed{
\lambda_i
=
\frac{w_i}{\sum_{j=1}^n w_j}
}
$$

so that

$$
\sum_i \lambda_i=1.
$$

If you want the **unnormalized version**, just use

$$
\boxed{
w_i
=
\frac{
\tan(\alpha_{i-1}/2)+\tan(\alpha_i/2)
}{
\sinh d_{\mathbb H}(p,q_i)
}
}
$$

directly.

Using disk coordinates only:

$$
\boxed{
w_i
=
\frac{
\tan(\alpha_{i-1}/2)+\tan(\alpha_i/2)
}{
\sinh\left(2\tanh^{-1}\left|\frac{q_i-p}{1-\overline p q_i}\right|\right)
}
}
$$

and since

$$
\sinh(2\tanh^{-1} \rho)=\frac{2\rho}{1-\rho^2},
$$

you can also write

$$
\boxed{
w_i
=
\left(\tan(\alpha_{i-1}/2)+\tan(\alpha_i/2)\right)
\frac{1-|u_i|^2}{2|u_i|}
}
$$

where

$$
u_i=\frac{q_i-p}{1-\overline p q_i}.
$$
