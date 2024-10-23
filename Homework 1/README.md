# MAE263F
Assignments for MAE263F FALL 2024

# Homework 1 Report

Contents
- `Homework1_LASTNAME.pdf`: Report addressing the questions for all three assignments, including plots and discussion.
- `Problem1.m`, `Problem2.m`, `Problem3.m` : MATLAB files for the respective problems.
- Helper Functions:
  - `gradEs.m`: Computes the gradient of the stretching energy.
  - `gradEb.m`: Computes the gradient of the bending energy.
  - `hessEs.m`: Computes the Hessian of the stretching energy.
  - `hessEb.m`: Computes the Hessian of the bending energy.

Instructions

Problem 1: Beam Dynamics - Implicit vs Explicit Simulation
- File: `Problem1.m`
- Description: This script simulates the dynamics of an elastic beam using both implicit and explicit numerical integration methods.
- How to Run:
  - Run the script in MATLAB.
  - You will be prompted to choose between the implicit or explicit simulation:
    - Enter `1` for the **implicit** method.
    - Enter `2` for the **explicit** method.
- Output The simulation will provide the beam's response based on the chosen method.

Problem 2: Mass-Spring System with Three Spheres
- File `Problem2.m`
- Description This script simulates a mass-spring system consisting of three spheres connected by an elastic beam. It models the vertical position and velocity of the middle sphere.
- How to Run
  - Run the script in MATLAB.
  - The output includes:
    - The shape of the beam over time.
    - The position and velocity of the second (middle) sphere.
- Additional Features
  - There are some commented sections within the code, allowing you to iterate through different time step (`dt`) and number of nodes (`N`) values for further analysis.

Problem 3: Comparison with Euler-Bernoulli Beam Theory
- File `Problem3.m`
- Description: This script simulates the deformation of an elastic beam subjected to a point load and compares the results with the theoretical prediction from Euler-Bernoulli beam theory.
- How to Run:
  - Run the script in MATLAB.
  - The output includes:
    - The shape of the beam over time.
    - The maximum vertical displacement (`ymax`) from the simulation.
    - A comparison between the simulated `ymax` and the predicted `ymax` from Euler-Bernoulli beam theory.
    - A plot illustrating the divergence between simulated and theoretical values for large deformations.
