from matplotlib import lines
import numpy as np 
from numpy import pi
import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

global i_anim

# Parse the results data for the pendulum results
# This parser assumes 4D phase space (2 angles and 2 rates)
# Inputs:
# - filename of the .csv file where the results are saved
# Output:
# It will return arrays of [ time , theta , phi , om_theta , om_phi ] with the same dimensionality (number of points)
# NOTE: In case a different system is used (more/less dimensions) the general parser bellow should be used instead
def parse_results_doublep( filename ):

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

# Parse the results data for the pendulum results
# This parser does not assume anything about the phase space and checks it dynamicslly
# Inputs:
# - filename of the .csv file where the results are saved
# Outputs:
# - time[ N ], x_res[ Nphase/2 ][ N ], v_res[ Nphase/2 ][ N ]
# NOTE: N is the time dimension, the x/v arrays have dimensionality 1/2 of the phase space
def parse_results_general( filename ):

    fp = open( filename , "r" )

    if fp.readable( ):
        data = csv.reader( fp )
        lst = [ ]
        for line in data:
            lst.append( line )
        ndata = len( lst ) - 1
        narr = ( len( lst[ 0 ] ) - 1 )/ 2
        print( narr )

        time = [ 0 ]*ndata
        x_res = np.zeros( ( narr , ndata ) )
        v_res = np.zeros( ( narr , ndata ) )

        for i in range( 0 , ndata ):
            time[ i ] = float( lst[ i + 1 ][ 0 ] )
            for j in range( 0 , narr ):
                x_res[ j ] = float( lst[ i + 1 ][ j ] )
                v_res[ j ] = float( lst[ i + 1 ][ j + narr ] )
    else:
        print( "Unreadable data, something's wrong with the file " + filename )
    
    return time, x_res, v_res

# Make a 2D plot of the angles and angular rates as functions of time
# Inputs:
# - time[ N ]: time array [sec] assumed
# - theta[ N ]: first angle array [rad] assumed
# - phi[ N ]: second angle array [rad] assumed
# - om_theta[ N ]: first angular rate [rad/s] assumed
# - om_phi[ N ]: second angular rate [rad/s] assumed
# - file_name: string which should contain "0" if no plot should be saved and the filename string if you want to save the figure in .pdf
# Output:
# - Show 2 plots (of the position and of the time) with the given assumed scales
# - Optionally save the figures based on the names provided
def plot_2D_time( time , theta , phi , om_theta , om_phi , file_name ):

    fig, ax = plt.subplots()
    ax.plot( time , theta , color = "green" , linestyle = "solid" , label = r"$\theta(t)$" )
    ax.plot( time , phi , color = "red" , linestyle = "solid" , label = r"$\varphi(t)$" )

    ax.set_xlabel( r"time [sec]" )
    ax.set_ylabel( r"Angles in [rad]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    if file_name != "0":
        fname = "Phi_Theta_Evolution_" + file_name + ".pdf"
        fig.savefig( fname , format = "pdf" )
    
    plt.show()  

    fig, ax = plt.subplots()
    ax.plot( time , om_theta , color = "green" , linestyle = "solid" , label = r"$\dot{\theta}(t)$" )
    ax.plot( time , om_phi , color = "red" , linestyle = "solid" , label = r"$\dot{\varphi}(t)$" )

    ax.set_xlabel( r"time [sec]" )
    ax.set_ylabel( r"Angular rates in [rad/s]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    if file_name != "0":
        fname = "dPhi_dTheta_Evolution_" + file_name + ".pdf"
        fig.savefig( fname , format = "pdf" )
    
    plt.show()  

# Make a 2D plot of the angles and angular rates in the phase space projections along the two x-v_x planes
# Inputs:
# - theta[ N ]: first angle array [rad] assumed
# - phi[ N ]: second angle array [rad] assumed
# - om_theta[ N ]: first angular rate [rad/s] assumed
# - om_phi[ N ]: second angular rate [rad/s] assumed
# - file_name: string which should contain "0" if no plot should be saved and the filename string if you want to save the figure in .pdf
# Output:
# - Show a plot (in the theta-om_theta and in the phi-om_phi cross-sections) with the given assumed scales
# - Optionally save the figures based on the names provided
def plot_2D_phase( theta , phi , om_theta , om_phi , file_name ):

    fig, ax = plt.subplots()
    ax.plot( theta , om_theta , color = "green" , linestyle = "solid" , label = r"$\omega_{\theta}(\theta)$" )
    ax.plot( phi , om_phi , color = "red" , linestyle = "solid" , label = r"$\omega_{\phi}(\phi)$" )

    ax.set_xlabel( r"Angles in [rad]" )
    ax.set_ylabel( r"Angular rates in [rad/s]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    if file_name != "0":
        fname = "Phase_Portrait_" + file_name + ".pdf"
        fig.savefig( fname , format = "pdf" )
    
    plt.show()  
