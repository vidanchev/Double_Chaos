#include <stdio.h>
#include <math.h>
#define PI 3.1415926536 /* I 8 sum ... and it was delicious */
#define Om 1.0 /* Angular frequency of the oscillator */

/* Some global variables which will be used for the integrators */
double DPc[ 7 ], /* Time-step coefficients {ci} from the Butcher Tableu */
       DPb[ 2 ][ 7 ], /* 4th and 5th order final weights {bi} and {b*i} from the Butcher Tableu */
       DPa[ 7 ][ 7 ], /* k computation weights {a_ij} -> Only the non-zero ones from the Butcher Tableu */
       DPec[ 7 ];  /* bi-b*i -> error coefficient constants (to avoid computation in each step) */

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
    - state_init[ dim_state ]: initial state for the integrator
    - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
    - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
    - header: a char of the header to start the file with (no need for \n sign) */
/* Outputs:
    - The results are written in a file as commas separated values (.csv)
    -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
    -- Reflect this in the header format! */
void RK4_Integrator( int Nstate , int Npoints , double* state_init , double* range_int , char* file_name , char* header ){

    int i, j; /* Iterators */
    double k_RK[ Nstate ][ 4 ], /* Runge-Kutta constants */
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

/* 4-5th order adaptive Dormand-Prince integrator - TO BE DONE */
/* Inputs:
    - Nstate: number of quantities in the state (phase space dimension)
    - state_init[ dim_state ]: initial state for the integrator
    - range_int[ 2 ]: initial and final values evolution parameter (initial and final time)
    - file_name: a char of the output filename where the results will be written - include .csv in this like "file.csv"
    - header: a char of the header to start the file with (no need for \n sign) */
/* Outputs:
    - The results are written in a file as commas separated values (.csv)
    -- The format is [ Time , State[ 0 ] , State[ 1 ] , ... State[ Nstate - 1 ] ]
    -- Reflect this in the header format! */
void DP45_Integrator( int Nstate , double* state_init , double* range_int , char* file_name , char* header ){

    int i, j; /* Iterators */
    double k_RK[ Nstate ][ 4 ], /* Runge-Kutta constants */
           state_now[ Nstate ], /* Variable where we keep the current state */
           int_state[ Nstate ], /* Variable where we keep intermediate state for RK steps */
           rhs_state[ Nstate ]; /* Right-Hand-Side of the state (derivatives) */ 
    double t_now, /* Current time value */
           dt; /* Time step */
    FILE *fp; /* File pointer to write the results */

    t_now = *( range_int ); /* Initialize time start */
    dt = ( *( range_int + 1 ) - *( range_int ) )/1e6; /* Initial "guess" for a good time step - will be changed from the integrator later */

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
    for( i = 0; i < 1; i++ ){
        
        printf( "No integration loop yet, you should write it :) \n" );

    }    

    fclose( fp ); /* Close the file in the end */

}

int main( ){

    int Npoints = 50;
    double state_init[ 2 ] = { 1.0 , 0.0 },
           range_int[ 2 ] = { 0.0 , 2.0*PI };
    char test_char[ ] = "Test_Results.csv";

    RK4_Integrator( 2 , Npoints , state_init , range_int , test_char , "T [time], X [pos], Y [vel]," );

    return 0;
}