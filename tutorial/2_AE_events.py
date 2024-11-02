"""
    - load AE lab data and compute fractal dimension
    - compute fractal dimension of arbitrary point cloud in 2D and 3D,
    compute correlation integral and fit linear portion in log-log space
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# add parent dir to path
p2parent = os.path.dirname( os.getcwd())
sys.path.append( p2parent)
from fractalAnalysis import fractalDim

#=================================0==================
#                set parameters
#====================================================
AE_file = "AE_PM_ST_152_002.txt" #"AE_PM_QZ_155_001.txt" #AE_PM_ST_152_002.txt
dPar = {
          # --- C(r) and D2----------
          'x_min'     : 0.1, #smallest distance for point density calculation and C(r)
           # sample length
           'length' : 30, # [cm]
            # max. fitting range
          'xmax'  : None, #3,
          #-------plotting-----------------------
          'testPlot'    : False,
          'plotFormat'  : 'png',
        }

#=================================1=========================================================
#                load data
#===========================================================================================
a_X, a_Y, a_Z = np.loadtxt( f"{p2parent}/data/{AE_file}", usecols=(1,2,3), skiprows=1,
                            dtype = float).T
N_AE = len( a_X)
if dPar['testPlot']:
    fig, ax = plt.subplots( 2,1, sharex=True)
    ax[0].set_title( f"N-AE={N_AE}")
    ax[0].plot( a_X, a_Y, 'ko', ms = 1)
    ax[0].set_xlabel('X (cm)'), ax[0].set_ylabel( 'Y (cm)')
    ax[0].set_xlim( 10, 40)
    ax[0].set_ylim( 0, 15)

    ax[1].plot( a_X, -a_Z, 'ko', ms = 1)
    ax[1].set_xlabel('X (cm)'), ax[0].set_ylabel( 'Z (cm)')
    ax[1].set_xlim(10, 40)
    ax[1].set_ylim( -10, 5)
    plt.show()

#=================================2=========================================================
#                    Correlation Integral
#===========================================================================================
d_C_int = fractalDim.C_r_Int( a_X, a_Y, a_Z,
                             x_min = dPar['x_min'],
                             x_max = dPar['length'],
                            logBinning = True, binFactor=1.1)
N_max = ( N_AE*(N_AE-1))*.5
# determine max. cut-off
if 'xmax' in dPar.keys() and dPar['xmax'] is not None:
    dRmax = { 'rmax' : dPar['xmax']}
else:
    dRmax = { 'rmax' : d_C_int['aL_ave'][d_C_int['aN_sum'] > N_max][0]}
print( 'Nmax:', N_max, dRmax)
sel_pl   = d_C_int['aL_ave'] < dRmax['rmax']
#=================================3=========================================================
#                        fractal dimension
#===========================================================================================
slope, interc, r, p, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][sel_pl]),
                                                  np.log10( d_C_int['aCorr'][sel_pl]))
interc = 10**interc
print( 'D2' , slope, 'rmax', dRmax['rmax'])
aC_hat = interc*d_C_int['aL_ave']**slope

color  = 'k'
plt.figure(1)
ax = plt.subplot( 111)
ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
            ms = 5, mec = color, mfc = 'w')
ax.loglog( d_C_int['aL_ave'], aC_hat, '--', color = color)
# mark upper bound
ax.axvline( dRmax['rmax'])
#------labels and legends
#ax.legend( loc = 'lower right')
ax.set_xlabel( 'Distance (cm)')
ax.set_ylabel( 'Correlation Integral')
plt.savefig( f"{p2parent}/plots/CorrInt_{AE_file.split('.')[0]}.{dPar['plotFormat']}")
plt.show()





    
    