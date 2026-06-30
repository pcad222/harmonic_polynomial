# Harmonic Polynomial Fitting for 3D Magnetic Field Reconstruction

This repository implements a harmonic-polynomial (solid spherical harmonic) method for reconstructing three-dimensional magnetic fields from discrete field measurements.

The magnetic field is represented as a linear combination of harmonic-polynomial basis functions. The expansion coefficients, $G_{\ell,m}$, are determined using a least-squares fit to the measured magnetic-field components.

The basis functions are derived from harmonic scalar potentials that satisfy Laplace's equation,

$$
\nabla^2 \Phi = 0,
$$

where the magnetic field is obtained from

$$
\mathbf{B} = -\nabla\Phi.
$$

Consequently, the reconstructed magnetic field automatically satisfies the magnetostatic Maxwell equations in a source-free region,

$$
\nabla \cdot \mathbf{B} = 0,
$$

and

$$
\nabla \times \mathbf{B} = 0.
$$

---

## Features

- Construct harmonic-polynomial basis functions up to degree `L_MAX`
- Compute the spherical harmonic coefficients $G_{\ell,m}$ using least-squares fitting
- Reconstruct the magnetic field from the fitted coefficients
- Compare reconstructed and measured magnetic fields
- Compute residual magnetic fields

$$
\mathbf{B}_{\rm residual}
=
\mathbf{B}_{\rm measured}
-
\mathbf{B}_{\rm reconstructed}
$$

- Visualize residual distributions using histograms for $B_x$, $B_y$, and $B_z$

---

## Repository Structure

```text
harmonic_polynomial/
│
├── data/
│   └── README.md
│
├── notebooks/
│   └── harmonic_polynomial_analysis.ipynb
│
├── src/
│   └── harmonic_polynomial.py
│
└── README.md
```

---

## Input Data

The analysis expects a magnetic-field dataset with the columns

```text
x, y, z, bx, by, bz
```

where

| Column | Description | Unit |
|--------|-------------|------|
| `x`, `y`, `z` | Measurement coordinates | m |
| `bx`, `by`, `bz` | Magnetic-field components | nT |

---

## Output

The notebook computes

- harmonic-polynomial coefficients $G_{\ell,m}$
- reconstructed magnetic field
- residual magnetic field
- residual histograms for $B_x$, $B_y$, and $B_z$

These results are used to evaluate the accuracy of the harmonic-polynomial reconstruction.

---

## Theory

The harmonic-polynomial expansion provides a physically constrained representation of the magnetic field by enforcing

- Laplace's equation
- $\nabla\cdot\mathbf{B}=0$
- $\nabla\times\mathbf{B}=0$

making it well suited for reconstructing magnetic fields in source-free regions.

---

## Dependencies

- NumPy
- Pandas
- SymPy
- Matplotlib
