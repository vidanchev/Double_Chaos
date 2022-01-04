# A script for running some tests before properly adding them to the different .py files

from Visualizations import parse_results_doublep
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

def animate( i ):

    x = theta[ i ]
    y = om_theta[ i ]
    x1.append( x )
    y1.append( y )

    x = phi[ i ]
    y = om_phi[ i ]
    x2.append( x )
    y2.append( y )

    xlist = [ x1 , x2 , [ theta[ i ] , phi[ i ] ] ]
    ylist = [ y1 , y2 , [ om_theta[ i ] , om_phi[ i ] ] ]

    #for index in range(0,1):
    for lnum,line in enumerate( lines ):
        line.set_data( xlist[ lnum ] , ylist[ lnum ] ) # set data for each line separately. 

    return lines

def init( ):
    for line in lines:
        line.set_data( [ ] , [ ] )
    return lines

# Load the file and parse the data
out_file = "Test_Results.csv" 
time, theta, phi, om_theta, om_phi = parse_results_doublep( out_file )

fig = plt.figure( )
ax1 = plt.axes( xlim = ( min( [ min( theta ) , min( phi ) ] ) , max( [ max( theta ) , max( phi ) ] ) ) , 
                ylim = ( min( [ min( om_theta ) , min( om_phi ) ] ) , max( [ max( om_theta ) , max( om_phi ) ] ) ) )
line, = ax1.plot( [ ] , [ ] , lw = 2 )
plt.xlabel( "Angle [rad]" )
plt.ylabel( "Angular Rates [rad/s]" )

plotlays, plotcols, plotleg = [ 3 ], [ "green" , "red" , "black" ], [ r"$\theta$" , r"$\varphi$" , "" ]
lines = []

for index in range( 3 ):
    lobj = ax1.plot( [ ] , [ ] , lw = 1 , color = plotcols[ index ] , label = plotleg[ index ] )[ 0 ]
    lines.append( lobj )

x1, y1 = [], []
x2, y2 = [], []

frame_num = len( time )

#x_arr = np.linspace( 0 , 2.0*np.pi , 100 )
#y_arr = np.sin( x_arr )
#z_arr = np.cos( x_arr )

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation( fig , animate , init_func = init ,
                               frames = frame_num , interval = 1 ) #, blit=True)

plt.show()

'''
from Visualizations import parse_results_doublep
from matplotlib import animation
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from numpy import pi

def init( ):
    ax.set_xlim( 0 , 2*np.pi )
    ax.set_ylim( - 1 , 1 )
    return ln,

def update( frame ):
    xdata.append( xarr[ frame ] )
    ydata.append( yarr[ frame ] )
    ln.set_data( xdata , ydata )
    return ln,


if __name__ == "__main__":

    out_file = b"Test_Results.csv" # Filename for the output - must be binary

    time, theta, phi, om_theta, om_phi = parse_results_doublep( out_file )
    
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln, = plt.plot([], [], 'r')

    xarr = np.linspace( 0 , 2*pi , 128 )
    yarr = np.sin( xarr )

    i_fr = [0]*128
    for i in range( 0 , 128 ):
        i_fr[ i ] = i

    ani = FuncAnimation( fig , update , frames = i_fr , init_func = init , blit = True , interval = 50 )

    #ani.save( "animation.gif" )

    plt.show()

'''

    #plot_2D_time( time , theta , phi , om_theta , om_phi , "0" )

    #plot_2D_phase( theta , phi , om_theta , om_phi , "0" )
    
    # A test for the harmonic oscillator call from the library as well -> it must result into cosine solution -> X(t) = \cos{t}
    
    #RK4_Integrator( 100 , [ 1.0 , 0.0 ] , [ 0.0 , 2.0*pi ] , b"Test_HO.csv" , b"T [time], X [pos], V [vel]" )