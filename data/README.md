## Data

This project uses **synthetic magnetic field data generated from simulations**. No experimental or measured data are included in this repository.

The input dataset is expected to be provided as a Pandas DataFrame with the following columns:

| Column | Description | Units |
|--------|-------------|-------|
| `x` | x-coordinate | m |
| `y` | y-coordinate | m |
| `z` | z-coordinate | m |
| `bx` | Magnetic field component along x | nT |
| `by` | Magnetic field component along y | nT |
| `bz` | Magnetic field component along z | nT |

Each row represents a single spatial measurement point.



The reconstruction algorithms assume that these six columns are present. Additional columns may be included but are ignored unless explicitly used.