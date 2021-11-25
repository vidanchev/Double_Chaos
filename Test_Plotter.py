import numpy as np 
import csv
import matplotlib.pyplot as plt
import os, sys

currentdir = os.path.dirname( os.path.realpath( __file__ ) )
parentdir = os.path.dirname( currentdir )
sys.path.append( parentdir )

# Parse the Test File data for the harmonic oscillator
def parse_data_HO( filename ):

    fp = open( filename , "r" )

    if fp.readable( ):
        data = csv.reader( fp )
        lst = [ ]
        for line in data:
            lst.append( line )
        ndata = len( lst ) - 1

        time = [ 0 ]*ndata
        pos = [ 0 ]*ndata 
        vel = [ 0 ]*ndata

        for i in range( 0 , ndata ):
            time[ i ] = float( lst[ i + 1 ][ 0 ] )
            pos[ i ] = float( lst[ i + 1 ][ 1 ] )
            vel[ i ] = float( lst[ i + 1 ][ 2 ] )
    else:
        print( "Unreadable data, something's wrong with the file " + filename )
    
    return time, pos, vel

# Parse the results data for the pendulum results
def parse_data_DP( filename ):

    fp = open( filename , "r" )

    if fp.readable( ):
        data = csv.reader( fp )
        lst = [ ]
        for line in data:
            lst.append( line )
        ndata = len( lst ) - 1

        time = [ 0 ]*ndata
        theta = [ 0 ]*ndata 
        phi = [ 0 ]*ndata
        om_theta = [ 0 ]*ndata
        om_phi = [ 0 ]*ndata

        for i in range( 0 , ndata ):
            time[ i ] = float( lst[ i + 1 ][ 0 ] )
            theta[ i ] = float( lst[ i + 1 ][ 1 ] )
            phi[ i ] = float( lst[ i + 1 ][ 2 ] )
            om_theta[ i ] = float( lst[ i + 1 ][ 3 ] )
            om_phi[ i ] = float( lst[ i + 1 ][ 4 ] )
    else:
        print( "Unreadable data, something's wrong with the file " + filename )
    
    return time, theta, phi, om_theta, om_phi

if __name__ == "__main__":

    file_name = "RK_C_Library/Test_Results.csv"

    time, theta, phi, om_theta, om_phi = parse_data_DP( file_name )

    # Testing against real solution for harmonic oscillator to make sure the integrator is OK
    om = 1.0
    t_real = np.linspace( 0 , 2.0*np.pi , 100 )
    x_real = np.cos( om*t_real )
    y_real = np.cos( 2.0*om*t_real ) 

    # First plot - the function
    fig, ax = plt.subplots()
    ax.plot( t_real , x_real , color = "blue" , linestyle = "solid" , label = r"Real Function at $\omega$" ) # Real function
    ax.scatter( time , theta , color = "orange" , label = r"Integrated Function at $\omega$" ) # Integrated Function
    ax.plot( t_real , y_real , color = "green" , linestyle = "solid" , label = r"Real Function at $2\omega$" ) # Real function
    ax.scatter( time , phi , color = "red" , label = r"Integrated Function at $2\omega$" ) # Integrated Function

    # Set labels, grid and ranges
    #ax.set( xlabel = "x" , ylabel = "y" , title = r"Stitched function $y(x)$" )
    #ax.grid()
    #plt.xlim( -10.0 , 10.0 )
    #plt.ylim( -2.0 , 2.0 )

    ax.legend( loc = "upper right" )

    #fig.savefig( "Function_Result.pdf" )
    plt.show()    

    err_th = [ 0 ]*len( time )
    err_phi = [ 0 ]*len( time )

    for i in range( 0 , len( time ) ):
        err_th[ i ] = abs( ( np.cos( om*time[ i ] ) - theta[ i ] ) )
        err_phi[ i ] = abs( ( np.cos( 2.0*om*time[ i ] ) - phi[ i ] ) )
    # Second plot - the error
    fig, ax = plt.subplots()
    ax.plot( time , err_th , color = "blue" , linestyle = "solid" , label = r"Absolute Error at $\omega$" ) 
    ax.plot( time , err_phi , color = "green" , linestyle = "solid" , label = r"Absolute Error at $2\omega$" )

    # Set labels, grid and ranges
    #ax.set( xlabel = "x" , ylabel = "y" , title = r"Stitched function $y(x)$" )
    #ax.grid()
    #plt.xlim( -10.0 , 10.0 )
    plt.ylim( - 5e-11 , 3e-10 )

    ax.legend( loc = "upper right" )

    #fig.savefig( "Function_Result.pdf" )
    plt.show() 