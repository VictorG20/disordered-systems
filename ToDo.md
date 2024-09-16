# To Do

How to reestructure the repository?

## 1 Belief Propagation

- What defines an Ising chain?
  - A number of nodes $N$, an external magnetic field $h$, and a coupling constant $J$.
- What defines a Belief Propagation problem on an Ising chain?
  - An Ising chain, an inverse thermal energy $\beta$, and inital probabilities for the boundaries, depending on what does one want to calculate.
  - We need to be able to compute left and right messages, each of these requires an initial probability or initial cavity precision.
  - In order to estimate the marginals we require probabilities from both sides.
- Create a plot that shows how, for increasing $N$, we get closer to the mean magnetization curve in the thermodynamic limit.
- Check documentation in html.
- Write tests

## 2.1 Direct Diagonalization
- Document all functions.

## 2.2 Cavity Method

- Multiply vector by sparse matrix?
- Sample Gaussian Wigner Matrix.
- Create a new group “cavity_method” within [exercise2.h5](data/exercise2.h5). Inside this group, a dataset must consist of:
        - The Gaussian Wigner Matrix as a sparse matrix, labeled by size and mean connectivity.
        - Eigenvalues of the matrix obtained by direct diagonalization, labeled by size and mean connectivity.
        - Datasets of cavity values, labeled by its epsilon value.
- `r, i = reim(z)`