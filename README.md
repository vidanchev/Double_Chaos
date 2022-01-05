# Double Chaos

## Introduction
The purpose of this repo is to prepare a small adaptive-step numerical solver in C which will be imported as a shared library in Python and applied to the double pendulum problem.
It's a project meant for students of differential equations and numerical methods to explore one of the simplest and most classical "chaotic systems".
The reason for writing the numerical solver in C is to make it faster as compared to a pure Python implementation.
I envision reusing it for future projects as well.

## Code Structure
The code and description is distributed in the following subfolders:

- **Main_Code** contains the main Python file using the C shared library, individual scripts for the runs and contains all the plotting functions:
    - **main.py** is the main code where a run parameters are defined and the integration + plotting is called, it also contains the animation for making the actual pendulum visualization (not the static plots).
    - **RK_Driver.py** performs all the ctypes casting and calls the shared library from **RK_C_Library** described bellow, it is imported in any other Py code.
    - **Visualizations.py** parses the result files and holds different visualizations (2D and 3D animations)
    - **Test_Environment.py** is just a script used to test some functionalities before properly structuring the Py files
- **Physics_Description** contains a LaTeX file which will be used to describe the physics of the problem and later contain some plots and results.
- **RK_C_Library** contains a C file and header file with adaptive step Runge-Kutta (Dormand Prince) implementation for the double pendulum problem: 
    - **RK_Library.c** contains the RK library which will be used for integration of the dynamical equations. Eventually this will be closed as a standalone library. Currently it contains a Dormand-Prince O(4-5) intrinsic adaptive method but a RK(4) was also used for verification purposes. 
- **Test_Plotter.py** contains some parsers and plotters for experimental files which will be used while verifying the RK library.

**NOTE:** Currently there are several testing .csv files which compare results with harmonic oscillator and between the original C code nad the now shared library, these will be deleted later but are used for verification purposes.

## Recompiling the RK library
Note that Python is used only for plotting and setting the initial conditions, the actual numerical integration happens in a shared library, contained in the **RK_C_Library** folder.
If you make any modifications to the library you will have to recompile it and reflect any differences in the Python driver to make sure that everything works.
This can be done either in a 2-step compillation (passing through binary file) or in 1-step
compillation (shown with gcc):
- Two-step compillation: cfile to binary and binary to shared library:
  
        gcc -c -fPIC <cfile.c> -o <binary.o> 
        gcc <binary.o> -shared -o <libname.so>

- One-step compillation:
        
        gcc -shared -o <libname.so> -fPIC <cfile.c>
Where <cfile.c> = **RK_Library.c** and <libname.so> should be the name which you use to import the shared library in the Python driver. 

**NOTE:** You may also have to recompile the library in case you are running on a different system. I am using Mac so the extension is **.so**, which is also valid for Linux, under Windows that would be a **.lib** file.

