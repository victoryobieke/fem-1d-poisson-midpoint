# fem-1d-poisson-midpoint
# Finite Element Method for 1D Poisson Problem

This MATLAB project solves the boundary value problem:

\[
- u''(x) = f(x), \quad x \in (0, 1), \quad u(0) = u(1) = 0
\]

using the **Finite Element Method (FEM)** with **piecewise linear basis functions** and **midpoint quadrature** for integration.

---

## 🔧 Features

- Implements FEM for a 1D second-order elliptic PDE.
- Uses midpoint rule for assembling the load vector.
- Computes the following error norms:
  - \( L^2 \) norm
  - \( H^1 \) semi-norm (energy norm)
  - Maximum nodal error
- Demonstrates convergence rates with decreasing mesh size.
- Visualizes:
  - Error convergence on log-log plots
  - FEM vs exact solution comparison

---

## 📁 Files

- `fem_1d_poisson.m`: Main script for assembling, solving, and evaluating the FEM approximation.

---

## 🧮 Problem Setup

- **Exact solution:** \( u(x) = \sin(\pi x) \)
- **Source term:** \( f(x) = \pi^2 \sin(\pi x) \)
- **Boundary conditions:** Dirichlet, \( u(0) = u(1) = 0 \)

---

## 📈 Output

- Table showing error values and convergence orders
- Plot of FEM and exact solution on the finest mesh
- Log-log plots of error norms with reference slopes

---

## 📦 Requirements

- MATLAB R2018b or later (no additional toolboxes required)

---

## 📜 License

This project is released under the MIT License.

---

## 🙋‍♂️ Author

Victory Obieke
Oregon State University
