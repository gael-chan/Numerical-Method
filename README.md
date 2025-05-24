# Numerical-Method
**Abstract**  
This project investigates numerical methods for pricing arithmetic Asian options and up-and-out barrier options. Analytical and Monte Carlo techniques are used to compute option prices and first-order sensitivities (deltas), including control variates and antithetic variates to enhance computational efficiency. MATLAB scripts implement exact formulas, Monte Carlo simulations, and convergence analysis. The work provides insights into variance reduction, estimator stability, and model convergence towards continuous-monitoring limits.

---

## 🧪 Coursework Scope

This repository contains the implementation and analysis for two major exercises:

### 📌 Exercise 1: Asian Option Pricing and Sensitivity Analysis
- **Exact Pricing** of:
  - Arithmetic Asian Option (Monte Carlo)
  - Lower Bound (Curran, 1994)
  - Geometric Asian Option (Kemna and Vorst, 1990)
- **Sensitivity Estimators (Delta)** using:
  - Likelihood Ratio method
  - Pathwise method (where applicable)
- **Control Variate Techniques**:
  - Lower Bound as control
  - Geometric option as control
- **Efficiency Evaluation**:
  - Standard error, confidence intervals, and runtime
  - Efficiency ratio:  
    $begin:math:display$
    E(K, n) = \\frac{t_a \\cdot \\sigma_a^2}{t_b \\cdot \\sigma_b^2}
    $end:math:display$

### 📌 Exercise 2: Barrier Option Simulation and Convergence
- **Up-and-Out Barrier Call Option**:
  - Monte Carlo with and without antithetic variates
  - Convergence study as number of monitoring dates increases
- **Benchmarking** against closed-form solution for continuous monitoring (Shreve, 2004)

---

## 📁 Repository Structure

```bash
numerical-methods-project/
├── Ex1_Q3_closed_formula.m                  # Exact formulas for lower bound and geometric option
├── Ex1_Q4_Crude_Simulation_Price.m          # Crude Monte Carlo pricing of Asian options
├── Ex1_Q4_Crude_Simulation_Sensitvity.m     # Crude simulation of delta
├── Ex1_Q4_price_control_variate.m           # Control variate pricing
├── Ex1_Q4_first_sensitivity_control_variate.m # Control variate delta estimation
├── Ex1_Q5_Convergence_Analysis.mlx          # Discrete vs Continuous LB convergence analysis
├── SMM313_CW_Report_Group4.pdf              # Final project report (with results, discussion)
├── Group coursework NM.pdf                  # Coursework requirements
├── README.md                                # Project overview (this file)
```

---

## 🔍 Key Insights

- Control variates significantly reduce variance:
  - **Lower Bound CV** reduces price estimation error by ~27x
  - **Geometric CV** achieves ~2.6x error reduction
- High correlation (ρ > 0.999) between CV and target estimator enables dramatic efficiency gain
- Discrete LB converges to Continuous LB with **O(1/n)** rate
- Antithetic variates reduce standard error by over 60%, even for path-dependent payoffs like barrier options



## 🧑‍💻 Authors

- Gael Chan  
- Huiqi Li  
- Nathaniel Stephenson  
- Wincy So  
*(Bayes Business School – MSc Quantitative Finance)*

---
