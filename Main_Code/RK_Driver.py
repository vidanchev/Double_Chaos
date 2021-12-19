# This driver is used to simplify the calls of the C shared library
# Its functions act as an interface for ctypes conversion:
# Simply import it into any other Python code (from the respective directory) and run the functions as described!

from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

# Import the shared RK library
lib_RK = CDLL( "../RK_C_Library/RK_Library.so" )

# IMPORTANT: Must be done before running any integration !!!
# Initialize all the RK coefficients for the Dormand-Prince method !!!
lib_RK.Set_RK_Coeff( )

# Prinout the Runge-Kutta constants for integration to check that they are set appropriately
def Check_RK_Coeff( ):

    lib_RK.Check_RK_Coeff( )

# Test interface to the C library from Py 
# Enter double x value to be allocated and check that it is true 
def Test_Clib_Interface( x_val ):

    lib_RK.Test_Interface( c_double( x_val ) )

# Set the pendulum coefficients for the RHS computation
# Inputs:
# - coeff_vals[ 5 ] double array contains the coefficients a_th to b_phi in sequence
# NOTE: This function must be called before starting an integration or all the constants will be defaulted to 1.0 
def Set_Pend_coeff( coeff_vals ):

    lib_RK.Set_Pend_coeff.restype = None
    lib_RK.Set_Pend_coeff.argtypes = [ ndpointer( c_double ) ]
    lib_RK.Set_Pend_coeff( np.array( coeff_vals ) )

# 4th order Runge-Kutta integrator for testing purposes
# Inputs:
# - Npoints: number of integration points (NOT INTERVALS)
# - state_init[ dim_state ]: initial state for the integrator
# - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
# - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
# - header: a char of the header to start the file with (no need for \n sign) */
# Outputs:
# - The results are written in a file as commas separated values (.csv)
# -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
# -- Reflect this in the header format!
# NOTE: Currently this function is linked to a RHS of a single Harmonic Oscillator with angular frequency \omega = 1.0
def RK4_Integrator( npoints , state_init , range_int , file_name , header ):

    nstate = len( state_init )

    lib_RK.RK4_Integrator.restype = None
    lib_RK.RK4_Integrator.argtypes = [ c_int , c_int , ndpointer( c_double ) , ndpointer( c_double ) , c_char_p , c_char_p ]
    lib_RK.RK4_Integrator( nstate , npoints , np.array( state_init ) , np.array( range_int ) , file_name , header )

# 4-5th order adaptive Dormand-Prince integrator 
# Inputs:
# - err_tol: error tolerance per step -> the adaptive step is modified to maintain this
# - state_init[ dim_state ]: initial state for the integrator
# - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
# - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
# - header: a char of the header to start the file with (no need for \n sign) 
# Outputs:
# - The results are written in a file as commas separated values (.csv)
# -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
# -- Reflect this in the header format!
def DP45_Integrator( err_tol , state_init , range_int , file_name , header ):

    nstate = len( state_init )

    lib_RK.DP45_Integrator.restype = None
    lib_RK.DP45_Integrator.argtypes = [ c_int , c_double , ndpointer( c_double ) , ndpointer( c_double ) , c_char_p , c_char_p ]
    lib_RK.DP45_Integrator( nstate , err_tol , np.array( state_init ) , np.array( range_int ) , file_name , header )