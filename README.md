# Computation of Klein-Gordon Eigenvalues
This repository contains a MATLAB routine that computes the eigenvalues of the Klein-Gordon equation in one space dimension. The algorithm is based on a recent research paper in the framework of the [Solvability Complexity Hierarchy](https://cordis.europa.eu/project/id/885904). The numerical method is based a Birman-Schwinger type decomposition of the Klein-Gordon operator, followed by a matrix approximation, which converges in operator norm. This process identifies the eigenvalues as the poles of a matrix-valued function $z\mapsto \\|(I-K(z))^{-1}\\|$ in the complex plane. The MATLAB routine returns a contour plot of $\\|(I-K(z))^{-1}\\|$, as well as the locations of its poles, computed via a gradient ascent method. The contents of this package are:
* `main.m`: Main script that computes and plots the eigenvalues for a potential given in `potential.m`;
* `potential.m`: Potential function;
* `build_K.m`: Computes the matrix approximation to the function $K(z)$;
* `simpson_integral.m`, `potential_integrals.m`: Technical functions that compute the matrix elements of the relevant operators via the Simpson's integral method;
* `GD.m`: Computes the poles of $\\|(I-K(z))^{-1}\\|$ via standard gradient ascent.

![a plot of the algorithm's output](https://github.com/frank-roesler/spectral_klein_gordon/blob/main/output.png)

Any comments or queries are welcome at https://frank-roesler.github.io/contact/
