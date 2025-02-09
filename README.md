# NJL Model Analysis of Dense Quark Matter

This project implements a numerical analysis of dense quark matter using a two-flavor NJL (Nambu–Jona-Lasinio) model approach. The code primarily calculates thermodynamic quantities such as the grand potential, quark number density, normalized baryon density, pressure, and energy density by solving a self-consistent gap equation for the effective quark mass \(M\) and an effective chemical potential \(\mu_{\text{tilde}}\).

## Features

- **Fermi Occupation Numbers:**  
  - `np(p, T, mu, M)`: Computes the Fermi-Dirac occupation number for quarks.
  - `np_bar(p, T, mu, M)`: Computes the Fermi-Dirac occupation number for antiquarks.

- **Thermodynamic Potential:**  
  - `OmegaM(Temp, μ_t, Mass)`: Computes the free (Fermi-gas) thermodynamic potential.
  - `Omega_temp(T, μ, M, μ_tilde)`: Computes the total thermodynamic potential including interaction terms.

- **Gap Equation and Self-Consistency:**  
  - `gap_eq(M; T, mu)`: Evaluates the gap equation (self-consistency condition) for the effective mass.
  - `find_M(Temp, μ)`: Numerically solves the gap equation to determine the effective mass \(M\).

- **System Solver for Equilibrium:**  
  - `solve_system(T, μ)`: Solves the coupled equations \(\frac{d\Omega}{dM}=0\) and \(\frac{d\Omega}{d\mu_{\text{tilde}}}=0\) using `NLsolve` to obtain the equilibrium values of \(M\) and \(\mu_{\text{tilde}}\).

- **Thermodynamic Quantities:**  
  - `Omega(T, μ)`: Computes the grand potential (relative to vacuum).
  - `num_dens(T, mu, M)`: Computes the total quark number density via 3D momentum integration.
  - `rhoB(T, μ)`: Computes the normalized baryon number density.
  - `pressure(T, μ)`: Computes the pressure.
  - `energy_dens(T, μ)`: Computes the energy density.

## Requirements

- **Julia 1.6 or later** is recommended.
- The project requires the following Julia packages:
  - `QuadGK` for numerical integration.
  - `NLsolve` for solving nonlinear systems.
- A file named `constants.jl` must be included in the project; it should define the global constants, for example:  
  `Nf_o`, `Nc_o`, `Lambda_o`, `G`, `m_o`, `hbarc`, and `rho0`.

## Installation and Setup

1. **Clone the Repository:**
   ```bash
   git clone <repository_url>
   cd <repository_directory>

