# Population Dynamics Algorithm

## Create an equilibrated population of cavity precisions

### Algorithm outline

1. Initialize a complex population of size $N_{p}$ for the cavity precisions: $\hat{P} = (\omega_{1}, \dots, \omega_{N_{p}})$. For the current purpose, we require that the cavity precisions start with $\text{Re}[\omega_{\ell}] > 0$.
1. For a maximum number of iterations:
    * Perform a *sweep* of the population.
    * Check if equilibrium has been reached. 
    * If equilibrium, stop. Else, continue.

### Sweep

A sweep consists of the following steps:
1. Pick $c - 1$ random elements $\left\{ \omega_{\ell} \right\}$ from $\hat{P}$.
1. Sample $E$ from $\rho(E)$.
1. Choose a random element of the population, $\omega_{k}$, and replace it with
    $$
        \omega_{k} = i (\lambda - i \epsilon - E) + \sum_{\ell = 1}^{c - 1} \frac{1}{\omega_{\ell}}.
    $$
A sweep is completed when every member of the population $\hat{P}$ has been updated once according to the previous steps.

### Checking equilibration

The aim is for the population $\hat{P}$ to reach a point where *the distribution is not changing significantly anymore*. This does **not** mean that the population is not changing, but that properties of the distribution, such as the moments, change less after every sweep.