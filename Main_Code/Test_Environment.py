# A script for running some tests before properly structuring all the Py files

from RK_Driver import Test_Clib_Interface, Check_RK_Coeff, Set_Pend_coeff, DP45_Integrator, RK4_Integrator
from numpy import pi

if __name__ == "__main__":

    Test_Clib_Interface( pi )

    err_tol = 1e-12 # Error Tolerance
    state_init = [ 1.0 , 1.0 , 0.0 , 0.0 ] # Initial state
    range_int = [ 0.0 , 10.0*pi ] # Range of integration
    pend_par = [ 0.006 , 0.0015 , 0.0045 , 0.441 , 0.147 ] # Test parameters corresponding to l_1 = l_2 = 0.3 m, m_1 = m_2 = 0.1 kg 
    out_file = b"Test_Results.csv" # Filename for the output - must be binary
    header = b"T [time], Theta [rad], Phi [rad], Om_Theta [rad/s], Om_Phi [rad/s]" # Header for the output file

    Set_Pend_coeff( pend_par )

    DP45_Integrator( err_tol , state_init , range_int , out_file , header )

    # A test for the harmonic oscillator call from the library as well -> it must result into cosine solution -> X(t) = \cos{t}
    
    RK4_Integrator( 100 , [ 1.0 , 0.0 ] , [ 0.0 , 2.0*pi ] , b"Test_HO.csv" , b"T [time], X [pos], V [vel]" )