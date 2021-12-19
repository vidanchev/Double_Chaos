/* ---------------------------------------------------------------------------------------------------- */
/* Library for numerical integration of ordinary differential equations, it includes adaptive-step 
Dormand Prince 4-5th order and fixed Runge-Kutta 4th order. The RHS functions control the dynamical problem,
the actual integrator is independent of them. Change the RHS functions to change the differential equation! */
/* ---------------------------------------------------------------------------------------------------- */
/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <vidanchev@uni-sofia.bg> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return. Victor Ivaylov Danchev
 * ----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#define PI 3.1415926536 /* I 8 sum ... and it was delicious */
#define Om 1.0 /* Angular frequency of the oscillator FOR TESTING PURPOSES */
#define Nloop_max 1e5 /* Maximum number of iterations for the Dormand-Prince loop regardless of step */
#define safe_fac 0.8 /* Safety factor when choosing a smaller step after rejected one */
#define p_gain ( - 0.14 ) /* Proportional gain for step decrease (in power of error ratio) */
#define p_loss ( - 0.2 ) /* Proportional "loss" for step increase (in power of error ratio) */
#define i_gain ( - 0.08 ) /* Integral gain for step decrease (in power of error ratio) */
#define step_mrat 8.0 /* Maximum ratio of the new step with respect to the previous one (increase) */

/* Some global variables which will be used for the integrators */
double DPc[ 7 ], /* Time-step coefficients {ci} from the Butcher Tableu */
       DPb[ 2 ][ 7 ], /* 4th and 5th order final weights {bi} and {b*i} from the Butcher Tableu */
       DPa[ 7 ][ 7 ], /* k computation weights {a_ij} -> Only the non-zero ones from the Butcher Tableu */
       DPec[ 7 ];  /* bi-b*i -> error coefficient constants (to avoid computation in each step) */

/* Some global variables which will be used for the pendulum characteristics */
double a_th = 1.0, /* a_{\theta} coefficient in the Lagrangian (term in front of \theta^2) */
       a_phi = 1.0, /* a_{\varphi} coefficient in the Lagrangian (term in front of \phi^2) */
       a_mix = 1.0, /* a_{\mix} coefficient in the Lagrangian (term in front of the mixed cosines) */
       b_th = 1.0, /* b_{\theta} coefficient in the Lagrangian  */
       b_phi = 1.0; /* b_{\varphi} coefficient in the Lagrangian  */

/* Populate the Runge-Kutta constants for integration
 - currently hardcoded for 4-5th order Dormand Prince adaptive step with embedded error estimation */
/* NOTE: Must be performed before any integrations with Dormand-Prince are performed!!! */
void Set_RK_Coeff( ){

    int i;

    /* Set Coefficients for Dormand-Prince Method (or another 4+5th order method) */
    /* Time-step coefficients {ci} */
    DPc[ 0 ] = 0.0;
    DPc[ 1 ] = 1.0/5.0;
    DPc[ 2 ] = 3.0/10.0;
    DPc[ 3 ] = 4.0/5.0;
    DPc[ 4 ] = 8.0/9.0;
    DPc[ 5 ] = 1.0;
    DPc[ 6 ] = 1.0;

    /* 4th order final weights {bi} */
    DPb[ 0 ][ 0 ] = 35.0/384.0;
    DPb[ 0 ][ 1 ] = 0.0;
    DPb[ 0 ][ 2 ] = 500.0/1113.0;
    DPb[ 0 ][ 3 ] = 125.0/192.0;
    DPb[ 0 ][ 4 ] = - 2187.0/6784.0;
    DPb[ 0 ][ 5 ] = 11.0/84.0;
    DPb[ 0 ][ 6 ] = 0.0;

    /* 5th order final weights {b*i} */
    DPb[ 1 ][ 0 ] = 5179.0/57600.0;
    DPb[ 1 ][ 1 ] = 0.0;
    DPb[ 1 ][ 2 ] = 7571.0/16695.0;
    DPb[ 1 ][ 3 ] = 393.0/640.0;
    DPb[ 1 ][ 4 ] = - 92097.0/339200.0;
    DPb[ 1 ][ 5 ] = 187.0/2100.0;
    DPb[ 1 ][ 6 ] = 1.0/40.0;

    /* Compute Error Coefficients */
    for( i = 0; i < 7; i++ ){
        DPec[ i ] = DPb[ 0 ][ i ] - DPb[ 1 ][ i ];
    }

    /* k computation weights -> Only the non-zero ones  */
    DPa[ 1 ][ 0 ] = 1.0/5.0;
    DPa[ 2 ][ 0 ] = 3.0/40.0;
    DPa[ 2 ][ 1 ] = 9.0/40.0;
    DPa[ 3 ][ 0 ] = 44.0/45.0;
    DPa[ 3 ][ 1 ] = - 56.0/15.0;
    DPa[ 3 ][ 2 ] = 32.0/9.0;
    DPa[ 4 ][ 0 ] = 19372.0/6561.0;
    DPa[ 4 ][ 1 ] = - 25360.0/2187.0;
    DPa[ 4 ][ 2 ] = 64448.0/6561.0;
    DPa[ 4 ][ 3 ] = - 212.0/729.0;
    DPa[ 5 ][ 0 ] = 9017.0/3168.0;
    DPa[ 5 ][ 1 ] = - 355.0/33.0;
    DPa[ 5 ][ 2 ] = 46732.0/5247.0;
    DPa[ 5 ][ 3 ] = 49.0/176.0;
    DPa[ 5 ][ 4 ] = - 5103.0/18656.0;
    DPa[ 6 ][ 0 ] = 35.0/384.0;
    DPa[ 6 ][ 1 ] = 0.0;
    DPa[ 6 ][ 2 ] = 500.0/1113.0;
    DPa[ 6 ][ 3 ] = 125.0/192.0;
    DPa[ 6 ][ 4 ] = - 2187.0/6784.0;
    DPa[ 6 ][ 5 ] = 11.0/84.0;

}

/* Prinout the Runge-Kutta constants for integration to check that they are set appropriately */
void Check_RK_Coeff( ){

    int i;

    printf( "DPc[ 7 ] coefficients are: " );
    for( i = 0; i < 7; i++ ){
        printf( "%.10e \n" , DPc[ i ] );
    }

    printf( "DPb[ 2 ][ 7 ] coefficients are ( DPb[ 0 ][ i ] , DPb[ 1 ][ i ] ): " );
    for( i = 0; i < 7; i++ ){
        printf( "%.10e, %.10e \n" , DPb[ 0 ][ i ] , DPb[ 1 ][ i ] );
    }    

    printf( "DPa[ 7 ][ 7 ] coefficients are ( DPa[ 0 ][ i ] , ... , DPa[ 6 ][ i ] ): " );
    for( i = 0; i < 7; i++ ){
        printf( "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e \n" , DPa[ 0 ][ i ] , DPa[ 1 ][ i ] , DPa[ 2 ][ i ] , DPa[ 3 ][ i ] , DPa[ 4 ][ i ] , DPa[ 5 ][ i ] , DPa[ 6 ][ i ] );
    } 

    printf( "DPec[ 7 ] coefficients are: " );
    for( i = 0; i < 7; i++ ){
        printf( "%.10e \n" , DPec[ i ] );
    }  

}

/* Test interface to the C library from Py */
/* Enter x value to be allocated and check that it is true */
void Test_Interface( double x_val ){

  printf( "If you're seeing this, the C function is running! \n" );
  printf( "The provided x value is %.10e \n" , x_val );

}

/* Set the dynamic problem coefficients */
/* Inputs:
    - coeff_vals[ 5 ] double array contains the coefficients a_th to b_phi in sequence */
/* NOTE: This function must be called before starting an integration or all the constants will be defaulted to 1.0 */
void Set_Pend_coeff( double *coeff_vals ){

    a_th = *( coeff_vals );
    a_phi = *( coeff_vals + 1 );
    a_mix = *( coeff_vals + 2 );
    b_th = *( coeff_vals + 3 );
    b_phi = *( coeff_vals + 4 );

}

/* Simple max value finder */
/* Inputs:
    - x[ N ] double array in which to find the max value
    - N = size of the array */
/* Output: 
    - res is the largest element of the array value */
double MaxVal( double *x , int N ){

    int i;
    double res, xnow;

    res = *( x );

    for( i = 1; i < N; i++ ){
        xnow = *( x + i );
        if( xnow >= res ){
            res = xnow;
        }
    }

    return res;
}

/* Right-Hand-Side Function for a 1D Harmonic Oscillator with angular rate Om defined above */
/* State is assumed to be [ position , velocity ] in arbitrary units */
/* Inputs:
    - state[ 2 ]: the state */
/* Outputs:
    - deriv_state[ 2 ]: the derivative of the state */
void RHS_Function_HO( double* state , double* deriv_state ){

    *( deriv_state ) = *( state + 1 );
    *( deriv_state + 1 ) = - Om*Om*( *( state ) );

}

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
void RK4_Integrator( int Nstate , int Npoints , double* state_init , double* range_int , char* file_name , char* header ){

    int i, j; /* Iterators */
    double k_RK[ Nstate ][ 4 ], /* Runge-Kutta intermediate derivatives */
           state_now[ Nstate ], /* Variable where we keep the current state */
           int_state[ Nstate ], /* Variable where we keep intermediate state for RK steps */
           rhs_state[ Nstate ]; /* Right-Hand-Side of the state (derivatives) */ 
    double t_now, /* Current time value */
           dt; /* Time step */
    FILE *fp; /* File pointer to write the results */

    t_now = *( range_int ); /* Initialize time start */
    dt = ( *( range_int + 1 ) - *( range_int ) )/( ( double )Npoints - 1.0 ); /* Get the interval size from number of points and total interval */

    /* Assign the initial state */
    for( j = 0; j < Nstate; j++ ){
        state_now[ j ] = *( state_init + j );
    }

    /* Open the file, write header and initial data */
    fp = fopen( file_name , "w" ); 
    fprintf( fp , "%s \n" , header );
    fprintf( fp , "%.10e, " , t_now );
    for( j = 0; j < Nstate; j++ ){
        fprintf( fp , "%.10e, " , state_now[ j ] );
    }
    fprintf( fp , "\n" );

    /* Start the main integration loop */
    for( i = 0; i < Npoints - 1; i++ ){
        
        /* Call the RHS function in the current point -> x_i */
        RHS_Function_HO( state_now , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( j = 0; j < Nstate; j++ ){
            k_RK[ j ][ 0 ] = rhs_state[ j ]*dt;
            int_state[ j ] = state_now[ j ] + k_RK[ j ][ 0 ]/2.0; /* x_i + k_1/2 */
        }

        /* Call the RHS function in the first half-point -> x_i + k_1/2 */
        RHS_Function_HO( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( j = 0; j < Nstate; j++ ){
            k_RK[ j ][ 1 ] = rhs_state[ j ]*dt;
            int_state[ j ] = state_now[ j ] + k_RK[ j ][ 1 ]/2.0; /* x_i + k_2/2 */
        }

        /* Call the RHS function in the second half-point -> x_i + k_2/2 */
        RHS_Function_HO( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( j = 0; j < Nstate; j++ ){
            k_RK[ j ][ 2 ] = rhs_state[ j ]*dt;
            int_state[ j ] = state_now[ j ] + k_RK[ j ][ 2 ]; /* x_i + k_3 */
        }

        /* Call the RHS function in the next point -> x_i + k_3 */
        RHS_Function_HO( int_state , rhs_state );
        /* Assign RK constant for this point - no more intermediate steps */
        for( j = 0; j < Nstate; j++ ){
            k_RK[ j ][ 3 ] = rhs_state[ j ]*dt;
        }

        /* All RK constants for the step have been populated - we can find the next step */
        for( j = 0; j < Nstate; j++ ){
            state_now[ j ] += ( k_RK[ j ][ 0 ] + 2.0*k_RK[ j ][ 1 ] + 2.0*k_RK[ j ][ 2 ] + k_RK[ j ][ 3 ] )/6.0;
        }
        t_now += dt; /* Increment time */

        /* Write the new state in the output file */
        fprintf( fp , "%.10e, " , t_now );
        //fprintf( fp , "%.10e, %.10e, %.10e, \n" , t_now , state_now[ 0 ] , state_now[ 1 ] );
        for( j = 0; j < Nstate; j++ ){
            fprintf( fp , "%.10e, " , state_now[ j ] );
        }
        fprintf( fp , "\n" );

    }    

    fclose( fp ); /* Close the file in the end */

}


/* Mockup RHS value for two decoupled harmonic oscillators - it was used for the initial testing to validate the DP integrator */
/* NOTE: The second oscillator (at 2*\omega) has dampening by a coefficient 2.0*\beta defined in the function */
/* State is assumed to be [ \theta , \phi , \omega_theta , \omega_phi ] in arbitrary units of angle and angular velocity */
/* Inputs:
    - state[ 4 ]: the state as [ theta , phi , om_theta , om_phi ] */
/* Outputs:
    - deriv_state[ 4 ]: the derivative of the state in the same order */
/*
void RHS_Function( double* state , double* deriv_state ){

    

    double bet = 0.2;

    *( deriv_state ) = *( state + 2 );
    *( deriv_state + 1 ) = *( state + 3 ); 
    *( deriv_state + 2 ) = - Om*Om*( *( state ) );
    *( deriv_state + 3 ) = - 4.0*Om*Om*( *( state + 1 ) ) - 2.0*bet*( *( state + 3 ) );

}
*/

/* Right-Hand-Side Function for the double Pendulum with properties defined above */
/* State is assumed to be [ \theta , \phi , \omega_theta , \omega_phi ] in arbitrary units of angle and angular velocity */
/* Inputs:
    - state[ 4 ]: the state as [ theta , phi , om_theta , om_phi ] */
/* Outputs:
    - deriv_state[ 4 ]: the derivative of the state in the same order */
void RHS_Function( double* state , double* deriv_state ){

    double invA[ 2 ][ 2 ]; /* Inverse of the LHS matrix A */
    double rhs_vec[ 2 ]; /* RHS vector */
    double sTh = sin( *( state ) ), /* \sin{ \theta } */
           cTh = cos( *( state ) ), /* \cos{ \theta } */
           sDel = sin( *( state + 1 ) - *( state ) ), /* \sin{ \phi - \theta } */
           cDel = cos( *( state + 1 ) - *( state ) ), /* \cos{ \phi - \theta } */
           detA = ( 4.0*a_th*a_phi - a_mix*a_mix*cDel*cDel ); /* Determinant of the Left-Hand-Side (LHS) to be inverted */
    
    /* Check to see if the system is close to singular and warn in that case */
    /* NOTE: Analythically this should not be possible but check anyway in case of error in the code! */
    if( fabs( detA ) < 1e-12 ){
        printf( "WARNING: The system appears to be going singular at the state: \n" );
        printf( "theta = %.10e, phi = %.10e, om_theta = %.10e, om_phi = %.10e \n" , *( state ) , *( state + 1 ) , *( state + 2 ) , *( state + 3 ) );
        printf( "Assuming detA == 1 to continue, results after this point are WRONG!" );
        detA = 1.0;
    }

    /* Get the inverse matrix of the LHS */
    invA[ 0 ][ 0 ] = 2.0*a_phi/detA;
    invA[ 0 ][ 1 ] = - a_mix*cDel/detA;
    invA[ 1 ][ 0 ] = - a_mix*cDel/detA;
    invA[ 1 ][ 1 ] = 2.0*a_th/detA;

    rhs_vec[ 0 ] = - b_th*sTh + a_mix*sDel*( *( state + 3 ) )*( *( state + 3 ) );
    rhs_vec[ 1 ] = - b_phi*sin( *( state + 1 ) ) - a_mix*sDel*( *( state + 2 ) )*( *( state + 2 ) );

    /* Write out the derivatives */
    *( deriv_state ) = *( state + 2 );
    *( deriv_state + 1 ) = *( state + 3 ); 
    *( deriv_state + 2 ) = invA[ 0 ][ 0 ]*rhs_vec[ 0 ] + invA[ 0 ][ 1 ]*rhs_vec[ 1 ];
    *( deriv_state + 3 ) = invA[ 1 ][ 0 ]*rhs_vec[ 0 ] + invA[ 1 ][ 1 ]*rhs_vec[ 1 ];

}

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
void DP45_Integrator( int Nstate , double err_tol , double* state_init , double* range_int , char* file_name , char* header ){

    int rej, /* rej is a rejection counter -> ( 0 , 1 ) if the last step was rejected */
        i, j, k; /* Iterators */
    double k_DP[ Nstate ][ 7 ], /* Dormand-Prince intermediate derivatives */
           state_now[ Nstate ], /* Variable where we keep the current state */
           int_state[ Nstate ], /* Variable where we keep intermediate state for RK steps */
           rhs_state[ Nstate ], /* Right-Hand-Side of the state (derivatives) */ 
           err_est[ Nstate ]; /* Estimated error for each of the state quantities */

    double t_now, /* Current time value */
           t_mid, /* Intermediate time step during the DP integration */
           dt, /* Current time step */
           err_ratio, /* Error ratio of actual to desired - current step */
           err_ratiOld, /* Error ratio of actual to desired - old step */
           int_frac, /* This is used for interpolation in case we overshoot the end of the interval */
           tv1; /* Temporary variables which can be reused to hold some intermediate computations */
    FILE *fp; /* File pointer to write the results */

    t_now = *( range_int ); /* Initialize time start */
    dt = ( *( range_int + 1 ) - *( range_int ) )/1e6; /* Initial "guess" for a good time step is 1 millionth of the interval - will be modified from the integrator when it starts */

    /* Assign the initial state */
    for( j = 0; j < Nstate; j++ ){
        state_now[ j ] = *( state_init + j );
    }

    /* Open the file, write header and initial data */
    fp = fopen( file_name , "w" ); 
    fprintf( fp , "%s \n" , header );
    fprintf( fp , "%.10e, " , t_now );
    for( j = 0; j < Nstate; j++ ){
        fprintf( fp , "%.10e, " , state_now[ j ] );
    }
    fprintf( fp , "\n" );

    k = 0; /* Zero-out the loop counter */
    err_ratiOld = 1.0; /* Initialize the "old" error fraction */

    /* Start the main integration loop */
    while( ( t_now < *( range_int + 1 ) ) && ( k < Nloop_max ) ){

        /* Call the RHS function in the current point -> x_i */
        RHS_Function( state_now , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 0 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k2 computation */
            int_state[ i ] = state_now[ i ] + DPa[ 1 ][ 0 ]*k_DP[ i ][ 0 ];
        }

        /* Call the RHS function in the second point -> x_i + c2*dt */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 1 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k3 computation */
            int_state[ i ] = state_now[ i ]; 
            /* start with the existing state and add all the contributions */
            for( j = 0; j < 2; j++ ){
                int_state[ i ] += DPa[ 2 ][ j ]*k_DP[ i ][ j ];
            }
        }

        /* Call the RHS function in the third point -> x_i + c3*dt */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 2 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k4 computation */
            int_state[ i ] = state_now[ i ]; 
            /* start with the existing state and add all the contributions */
            for( j = 0; j < 3; j++ ){
                int_state[ i ] += DPa[ 3 ][ j ]*k_DP[ i ][ j ];
            }
        }

        /* Call the RHS function in the forth point -> x_i + c4*dt */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 3 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k5 computation */
            int_state[ i ] = state_now[ i ]; 
            /* start with the existing state and add all the contributions */
            for( j = 0; j < 4; j++ ){
                int_state[ i ] += DPa[ 4 ][ j ]*k_DP[ i ][ j ];
            }
        }

        /* Call the RHS function in the fifth point -> x_i + c5*dt */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 4 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k6 computation */
            int_state[ i ] = state_now[ i ]; 
            /* start with the existing state and add all the contributions */
            for( j = 0; j < 5; j++ ){
                int_state[ i ] += DPa[ 5 ][ j ]*k_DP[ i ][ j ];
            }
        }

        /* Call the RHS function in the sixth point -> x_i + c6*dt */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point and compute intermediate state */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 5 ] = rhs_state[ i ]*dt;
            /* This intermediate state will be used for k7 computation */
            int_state[ i ] = state_now[ i ]; 
            /* start with the existing state and add all the contributions */
            for( j = 0; j < 6; j++ ){
                int_state[ i ] += DPa[ 6 ][ j ]*k_DP[ i ][ j ];
            }
        }

        /* Call the RHS function in the seventh point -> x_i + c7*dt -> last one (t + dt) */
        RHS_Function( int_state , rhs_state );
        /* Assign RK constant for this point */
        for( i = 0; i < Nstate; i++ ){
            k_DP[ i ][ 6 ] = rhs_state[ i ]*dt;
        }

        /* We have all the k_DP at this point -- compute the error estimates for each state quantity */
        for( i = 0; i < Nstate; i++ ){
            err_est[ i ] = 0.0;
            for( j = 0; j < 7; j++ ){
                err_est[ i ] += DPec[ j ]*k_DP[ i ][ j ];
            }
            /* Take absolute value for each error */
            err_est[ i ] = fabs( err_est[ i ] );
        }

        /* Get the largest estimated error and find its ratio to the tolerance */
        err_ratio = MaxVal( err_est , Nstate )/err_tol;

        /* Choose whether to accept the step or not and how to pick the next step based on the err_ratio */
        /* NOTE: If the err_ratio is OK but we overshot the endpoint by more than err_tol, reject the step with new dt to end on it exactly! */
        if( err_ratio < 1.0 && ( t_now + dt - *( range_int + 1 ) < err_tol ) ){
            /* In this case the step is small enough -> we're within the error range and not yet at the end of the interval */

            /* In this case compute the next values and print them in the file */
            t_now += dt;

            /* NOTE: No need to check if we overshot since this is included in the outer if statement! */

            /* Update the states based on the DP coefficients of 4th order (final weights) */
            for( i = 0; i < Nstate; i++ ){
                for( j = 0; j < 7; j++ ){
                    state_now[ i ] += DPb[ 0 ][ j ]*k_DP[ i ][ j ];
                }
            }

            /* If last step was not rejected - increase the current step with a safety factor based on the integral controller */
            if( rej == 0 && ( err_ratio > 0.0 ) ){
                /* Check what the new step candidate is and keep it in tv1 */
                tv1 = safe_fac*dt*pow( err_ratio , p_gain )*pow( err_ratiOld , i_gain );
                /* If step is increased more than step_mrat, increase it by step_mrat */
                if( tv1/dt < step_mrat ){
                    dt = tv1;
                }
                else{
                    dt *= step_mrat;
                }
            }

            /* If the step was accepted -> set the rejection ratio to 0 */
            rej = 0;

            /* Write the new state in the output file */
            fprintf( fp , "%.10e, " , t_now );
            for( i = 0; i < Nstate; i++ ){
                fprintf( fp , "%.10e, " , state_now[ i ] );
            }
            fprintf( fp , "\n" );

        }
        else{
            /* In this case the step is too large or we overshot -> we must reduce it based on the estimate and threshold OR based on interval */
            if( t_now + dt - *( range_int + 1 ) > err_tol ){
                /* In this case we overshot the last step - set it to end exactly at the end of the interval */
                dt = ( *( range_int + 1 ) - t_now );
                rej = 0; /* If we came to this loop because of dt overshooting, we should not reject our new step */
            }
            else{
                /* In this case we're still integrating, reduce the step */
                dt = safe_fac*dt*pow( err_ratio , p_loss );
                rej = 1; /* Set rej to 1 in case the step was rejected */
            }

        }

        /* After a step has been completed - assign the new error ratio as old for next step */
        err_ratiOld = err_ratio; /* NOTE: This is used as an integral component in the step control */
        k += 1;

        /* In case we reached the maximum number of iterations - warn about it */
        if( k == Nloop_max ){
            printf( "----------------------------------------------------------\n" );
            printf( "----WARNING: The full integration was not carried out!----\n" );
            printf( "----------------------------------------------------------\n" );
            printf( "Stopped after %d iterations at t = %lf out of t_max = %lf \n" , k , t_now , *( range_int + 1 ) );

        }

        printf( "At k = %d, Time = %.10e, Error ratio is %.10e, dt = %.10e \n" , k , t_now , err_ratio , dt );


    }

    fclose( fp ); /* Close the file in the end */

}