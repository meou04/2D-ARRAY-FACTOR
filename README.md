# 2D Array Factor (AF) Simulation for Planar Antenna Arrays

## 📌 Project Overview
This repository contains a **MATLAB-based simulation** tool designed to calculate and visualize the **2D Array Factor (AF)** for planar (rectangular) antenna arrays. This project is essential for understanding phased array systems, beamsteering, and radiation pattern synthesis in 5G and satellite communications.

## 🚀 Key Features
* **N x M Planar Array Configuration:** Fully customizable number of elements and inter-element spacing ($d_x, d_y$).
* **Beamsteering Logic:** Calculation of progressive phase shifts ($\beta_x, \beta_y$) to steer the main lobe to any $(\theta_0, \phi_0)$ direction.
* **Side Lobe Control:** Implementation of amplitude tapering (Uniform, Binomial, or Taylor) to optimize the Side Lobe Level (SLL).
* **Advanced Visualization:** * 3D Spherical Radiation Patterns.
    * 2D Polar and Rectangular plots for E-plane and H-plane analysis.

## 📐 Mathematical Foundation
The Array Factor for a planar array in the $xy$-plane is computed using:

$$AF(\theta, \phi) = \sum_{m=1}^{M} \sum_{n=1}^{N} I_{mn} e^{j \left[ (m-1)(kd_x\sin\theta\cos\phi + \beta_x) + (n-1)(kd_y\sin\theta\sin\phi + \beta_y) \right]}$$

Where:
* $k = 2\pi/\lambda$ (Wavenumber).
* $I_{mn}$ is the excitation amplitude.
* $\beta_x, \beta_y$ are the phase shifts for beamsteering.

## 💻 Technical Stack
* **Language:** MATLAB
* **Required Toolboxes:** Antenna Toolbox (optional, if using built-in functions) or Signal Processing Toolbox.

## 📊 Sample Results
*(Tip: Upload your simulation screenshots to an 'images' folder and link them below)*
![3D Plot](images/your_simulation_result.png)
*Example: $8 \times 8$ Planar Array Simulation with $\lambda/2$ spacing.*

## 🛠 How to Run
1. Clone the repository:
   ```bash
   git clone [https://github.com/](https://github.com/)[Your-Username]/[Your-Repo-Name].git
