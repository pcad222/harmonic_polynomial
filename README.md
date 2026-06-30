# Harmonic Polynomial Fitting for 3D Magnetic Field Reconstruction

This repository implements a harmonic-polynomial method for reconstructing three-dimensional magnetic fields from discrete field measurements.

The magnetic field is represented as a linear combination of harmonic-polynomial basis functions. The best-fit coefficients, **Gₗₘ**, are computed using least-squares fitting.

The basis functions are derived from harmonic scalar potentials that satisfy Laplace's equation. Therefore, in a source-free region, the reconstructed magnetic field satisfies:

```text
div B = 0
curl B = 0
```

This makes the method physically constrained and suitable for magnetic-field reconstruction.

---

## Features

- Build harmonic-polynomial basis functions up to degree `L_MAX`
- Compute best-fit **Gₗₘ** coefficients using least-squares fitting
- Reconstruct the magnetic field from the fitted coefficients
- Compare the original and reconstructed magnetic fields
- Compute residuals using

```text
Residual = Truth − Predicted
```

- Plot residual histograms for `bx`, `by`, and `bz`

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

The input data should contain the following columns:

```text
x, y, z, bx, by, bz
```

where:

- `x`, `y`, `z` are measurement coordinates in meters
- `bx`, `by`, `bz` are magnetic-field components in nT

---

## Main Workflow

1. Load magnetic-field measurement data
2. Construct harmonic-polynomial basis functions
3. Fit the **Gₗₘ** coefficients using least squares
4. Reconstruct the magnetic field
5. Compute residuals between the truth and predicted fields
6. Plot residual distributions

---

## Output

The notebook computes:

- fitted **Gₗₘ** coefficients
- reconstructed magnetic-field components
- residual magnetic-field components
- residual histograms for `bx`, `by`, and `bz`

These outputs are used to evaluate the accuracy of the harmonic-polynomial reconstruction.

---

## Dependencies

- NumPy
- Pandas
- Matplotlib

---

## Notes

This repository is intended for research and analysis of magnetic-field reconstruction using harmonic-polynomial basis functions.
