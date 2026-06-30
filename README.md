# Harmonic Polynomial Fitting for 3D Magnetic Field Reconstruction

This repository implements a harmonic-polynomial (solid spherical harmonic) method for reconstructing three-dimensional magnetic fields from discrete field measurements.

The magnetic field is represented as a linear combination of harmonic-polynomial basis functions. The expansion coefficients, \(G_{\ell,m}\), are obtained using a least-squares fit to the measured magnetic-field components.

Because the basis functions are derived from solutions of Laplace's equation,

\[
\nabla^2 \Phi = 0,
\]

the reconstructed magnetic field

\[
\mathbf{B} = -\nabla\Phi
\]

automatically satisfies the magnetostatic Maxwell equations in a source-free region,

\[
\nabla\cdot\mathbf{B}=0,
\]

and

\[
\nabla\times\mathbf{B}=0.
\]

---

## Features

- Constructs harmonic-polynomial basis functions up to degree `L_MAX`
- Computes spherical harmonic coefficients \(G_{\ell,m}\) using least-squares fitting
- Reconstructs the magnetic field from the fitted coefficients
- Compares reconstructed and measured magnetic fields
- Computes residuals

\[
\mathrm{Residual}
=
\mathbf{B}_{\mathrm{measured}}
-
\mathbf{B}_{\mathrm{reconstructed}}
\]

- Produces residual histograms for \(B_x\), \(B_y\), and \(B_z\)

---

## Repository Structure

```text
harmonic_polynomial/
│
├── data/
│   └── README.md
│
├── src/
│   └── harmonic_polynomial.py
│
├── notebooks/
│   └── harmonic_polynomial_analysis.ipynb
│
└── README.md
```

---

## Data

The analysis expects magnetic-field measurements with the columns

```text
x, y, z, bx, by, bz
```

where

- `x`, `y`, `z` are measurement coordinates (m)
- `bx`, `by`, `bz` are magnetic-field components (nT)

---

## Output

The notebook computes

- fitted \(G_{\ell,m}\) coefficients
- reconstructed magnetic field
- residual distributions
- residual histograms

which are used to evaluate the quality of the harmonic-polynomial reconstruction.