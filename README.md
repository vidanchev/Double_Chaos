# Double_Chaos

## Introduction
The purpose of this repo is to prepare a small adaptive-step numerical solver in C which will be imported as a shared library in Python and applied to the double pendulum problem.
It's a project meant for students of differential equations and numerical methods to explore one of the simplest and most classical "chaotic systems".
The reason for writing the numerical solver in C is to make it faster as compared to a pure Python implementation.
I envision reusing it for future projects as well.

## Code Structure
The code and description is distributed in the following subfolders:

- **Main_Code** contains the Python driver for the C shared library, individual scripts for the runs and contains all the plotting functions, it also contains a few test scripts for experimentation: (**TO BE FILLED**)
- **Physics_Description** contains a LaTeX file which will be used to describe the physics of the problem and later contain some plots and results
- **RK_C_Library** contains a C file with adaptive step Runge-Kutta (Dormand Prince) implementation for the double pendulum problem: (**TO BE FILLED**)

## Compiling the RK library
**TO BE CONTINUED**