# Cavity Method

* For any given $\lambda$ we can compute $\rho(\lambda)$ according to
$$
    \rho(\lambda) = \lim_{\epsilon \to 0} \frac{1}{\pi N} \sum_{j = 1}^{N} \text{Im}\left[\frac{i}{\omega_{j}}\right] = \lim_{\epsilon \to 0} \frac{1}{\pi N} \sum_{j = 1}^{N} \left( \frac{\text{Re}(\omega_{j})}{\lvert \omega_{j} \rvert^{2}} \right)
$$

* A system of size $N$ has $N$ marginal precisions: $\omega_{1}, \dots, \omega_{N}$ which can be computed according to
$$
    \omega_{j} = i (\lambda - i \epsilon - E_{j}) + \sum_{k \in \partial j} \frac{1}{\omega_{k}^{(j)}}.
$$
Then, for a mean connectivity $c$, there are $N \cdot c$ cavity precisions:
$$
    \omega_{k}^{(j)} = i (\lambda - i \epsilon - E_{k}) + \sum_{\ell \in \partial k \backslash j} \frac{1}{\omega_{\ell}^{(k)}}.
$$

Starting from the equation for the cavity precisions, how to achieve the fixed point with a given accuracy?
1. Start all $\omega_{k}^{(j)}$ with random complex values.
1. 

## Notes

What takes the most time? Getting the cavity precisions up to a desired precision. Afterwards, computing the marginals and the spectral density is straightforward.