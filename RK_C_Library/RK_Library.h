/* ---------------------------------------------------------------------------------------------------- */
/* Header File for the library for numerical integration of ordinary differential equations RK_Library.c
Modify accordingly if you are modifying the RK_Library.c RHS functions to change the dynamical equations */
/* ---------------------------------------------------------------------------------------------------- */
/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <vidanchev@uni-sofia.bg> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return. Victor Ivaylov Danchev
 * ----------------------------------------------------------------------------
 */

/* Test interface to the C library from Py */
/* Enter x value to be allocated and check that it is true */
void Test_Interface( double x_val );

/* Populate the Runge-Kutta constants for integration
 - currently hardcoded for 4-5th order Dormand Prince adaptive step with embedded error estimation */
/* NOTE: Must be performed before any integrations with Dormand-Prince are performed!!! */
void Set_RK_Coeff( );

/* Prinout the Runge-Kutta constants for integration to check that they are set appropriately */
void Check_RK_Coeff( );

/* Set the dynamic problem coefficients */
/* Inputs:
    - coeff_vals[ 5 ] double array contains the coefficients a_th to b_phi in sequence */
/* NOTE: This function must be called before starting an integration or all the constants will be defaulted to 1.0 */
void Set_Pend_coeff( double *coeff_vals );

/* Simple max value finder */
/* Inputs:
    - x[ N ] double array in which to find the max value
    - N = size of the array */
/* Output: 
    - res is the largest element of the array value */
double MaxVal( double *x , int N );

/* Right-Hand-Side Function for a 1D Harmonic Oscillator with angular rate Om defined above */
/* State is assumed to be [ position , velocity ] in arbitrary units */
/* Inputs:
    - state[ 2 ]: the state */
/* Outputs:
    - deriv_state[ 2 ]: the derivative of the state */
void RHS_Function_HO( double* state , double* deriv_state );

/* 4th order Runge-Kutta integrator for testing purposes */
/* Inputs:
    - Nstate: number of quantities in the state (phase space dimension)
    - Npoints: number of integration points (NOT INTERVALS)
    - state_init[ Nstate ]: initial state for the integrator
    - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
    - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
    - header: a char of the header to start the file with (no need for \n sign) */
/* Outputs:
    - The results are written in a file as commas separated values (.csv)
    -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
    -- Reflect this in the header format! */
void RK4_Integrator( int Nstate , int Npoints , double* state_init , double* range_int , char* file_name , char* header );

/* Mockup RHS value for two decoupled harmonic oscillators - it was used for the initial testing to validate the DP integrator */
/* NOTE: The second oscillator (at 2*\omega) has dampening by a coefficient 2.0*\beta defined in the function */
/* State is assumed to be [ \theta , \phi , \omega_theta , \omega_phi ] in arbitrary units of angle and angular velocity */
/* Inputs:
    - state[ 4 ]: the state as [ theta , phi , om_theta , om_phi ] */
/* Outputs:
    - deriv_state[ 4 ]: the derivative of the state in the same order */
/*
void RHS_Function( double* state , double* deriv_state );
*/

/* Right-Hand-Side Function for the double Pendulum with properties defined above */
/* State is assumed to be [ \theta , \phi , \omega_theta , \omega_phi ] in arbitrary units of angle and angular velocity */
/* Inputs:
    - state[ 4 ]: the state as [ theta , phi , om_theta , om_phi ] */
/* Outputs:
    - deriv_state[ 4 ]: the derivative of the state in the same order */
void RHS_Function( double* state , double* deriv_state );

/* 4-5th order adaptive Dormand-Prince integrator */
/* Inputs:
    - Nstate: number of quantities in the state (phase space dimension)
    - err_tol: error tolerance per step -> the adaptive step is modified to maintain this
    - state_init[ Nstate ]: initial state for the integrator
    - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
    - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
    - header: a char of the header to start the file with (no need for \n sign) */
/* Outputs:
    - The results are written in a file as commas separated values (.csv)
    -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
    -- Reflect this in the header format! */
void DP45_Integrator( int Nstate , double err_tol , double* state_init , double* range_int , char* file_name , char* header );
