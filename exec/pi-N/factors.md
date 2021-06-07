# Distribution of kinematical and other factors in the code

The process $\pi N \to Ne^+e^-$ is calculated using a factorization of the process: creation of virtual photon + decay of virtual photon. The differential cross section can be written as

$$
    d\sigma = \frac{1}{32(2\pi)^4|\mathbf{q}|s} 
    dM\,|\mathbf{k}|\,|\mathbf{k}_1| \sum_{\lambda,\lambda'}
    \left[d(\cos\theta_k) \rho^\text{cre}(\lambda',\lambda)\right]
    \left[d\Omega_{k_1} \rho^\text{dec}(\lambda,\lambda')\right].
$$

This we write as

$$
    d\sigma = 
    dM\frac{M}{\pi}\sum_{\lambda,\lambda'}
    \underbrace{ \frac{1}{32\pi s}\frac{|\mathbf{k}|}{|\mathbf{q}|}
                 d(\cos\theta_\mathbf{k})\, \rho^\text{cre}(\lambda',\lambda)
               }_{\gamma^* \text{creation}}
               \times
    \underbrace{ \frac{1}{16\pi^2}\frac{|\mathbf{k}_1|}{M}
                d\Omega_{k_1}\, \rho^\text{dec}(\lambda,\lambda')
               }_{\gamma^* \text{decay}},
$$
where the "$\gamma^*$ creation" part contains all factors appearing in the differential cross section of the process $\pi N\to N\gamma^*$, while the "$\gamma$ decay" part contains only the factors coming from the phase-space integral, including a factor of $(2\pi)^4$.

