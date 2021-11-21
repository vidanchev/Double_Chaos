import numpy as np 
import csv
import matplotlib.pyplot as plt
import os, sys

currentdir = os.path.dirname( os.path.realpath( __file__ ) )
parentdir = os.path.dirname( currentdir )
sys.path.append( parentdir )

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


if __name__ == "__main__":

    file_name = "RK_C_Library/Test_Results.csv"

    time, pos, vel = parse_data_HO( file_name )

    # Testing against real solution for harmonic oscillator to make sure the integrator is OK
    om = 1.0
    t_real = np.linspace( 0 , 2.0*np.pi , 100 )
    x_real = np.cos( om*t_real )

    # First plot - the function
    fig, ax = plt.subplots()
    ax.plot( t_real , x_real , color = "blue" , linestyle = "solid" , label = r"Real Function" ) # Real function
    ax.scatter( time , pos , color = "red" , label = r"Integrated Function" ) # Integrated Function

    # Set labels, grid and ranges
    #ax.set( xlabel = "x" , ylabel = "y" , title = r"Stitched function $y(x)$" )
    #ax.grid()
    #plt.xlim( -10.0 , 10.0 )
    #plt.ylim( -2.0 , 2.0 )

    ax.legend( loc = "upper right" )

    #fig.savefig( "Function_Result.pdf" )
    plt.show()    