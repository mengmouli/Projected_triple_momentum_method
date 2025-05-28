# Projected Triple Momentum Method in MATLAB

This repository provides MATLAB implementations of the Projected Triple Momentum Method proposed in the following paper:  
https://arxiv.org/abs/2503.13965

---

## **Overview**

The code implements and compares the following algorithms:

- **Projected Gradient Descent**
- **Projected Triple Momentum Method with Euclidean Projection (Algorithm 1)**
- **Projected Triple Momentum Method with Weighted-Norm Projection (Algorithm 2)**

The algorithms in this repository are tested on the following constrained quadratic optimization problem:

Minimize &nbsp;&nbsp;&nbsp;&nbsp;**f(x) = 0.5 * x' * F * x + p' * x**
&nbsp;&nbsp; where  &nbsp;&nbsp;**F = [100, -1; -1, 1]**   &nbsp;&nbsp;&nbsp;&nbsp;**p = [1; 10]**

Subject to:
&nbsp;&nbsp;&nbsp;&nbsp;**x' * Q * x <= 1**   &nbsp;&nbsp; with &nbsp;&nbsp;**Q = [1, 0; 0, 2]**

---

## **Tested Environment**

- **MATLAB Version**: R2024a Update 6 (`24.1.0.2689473`)
- **Toolboxes**: Symbolic Math Toolbox
- **Operating System**: Windows 11 Pro (`Build 26100`)
- **Java**: 1.8.0_202 (Oracle HotSpot 64-bit)
- **YALMIP**: 20230622
- **MOSEK**: 10.2.11


---

## **File Descriptions**

| File                              | Description                                                                                   |
| --------------------------------- | --------------------------------------------------------------------------------------------- |
| `Algorithm_1.m`                  | Runs the first projected triple momentum algorithm (Algorithm 1) and visualizes convergence.  |
| `Algorithm_2.m`                  | Runs Algorithm 2 using a Lyapunov matrix `P` constructed via IQC theory.                          |
| `Projected_Gradient_Algorithm.m` | Runs projected gradient descent to highlight solver precision limitations.                    |
| `Comparison_projg_alg_1_2.m`     | Compares the convergence behavior of all three algorithms on the constrained optimization problem.             |


## **Supporting Functions**

The following files are utility scripts or subfunctions required by the main algorithms. The main scripts automatically call them.

| File                                          | Purpose                                                                                      |
|-----------------------------------------------|----------------------------------------------------------------------------------------------|
| `analyze_quadratic_function.m`                | Computes the gradient, strong convexity constant `m`, and Lipschitz constant `L` for a given quadratic function. |
| `build_projection.m`                          | Builds a projection operator from a user-defined constraint function.                        |
| `project_yalmip.m`                            | Solves the convex projection problem using YALMIP and MOSEK.                                |
| `proj_gradient.m`                             | Runs the projected gradient descent algorithm given system matrices and projection.          |
| `lyapunov_matrix_IQC_triple_momentum_unconstrained.m` | Constructs a Lyapunov function matrix using IQC theory for analyzing Algorithm 2.                   |
| `P.mat`                                       | Precomputed Lyapunov matrix (used in one of the steps of Algorithm 2).                         |
