"""
    - create a 1D cantor set or
      2D sierpinski triangle and
      compute fractal dimension
    - note that this script has been
      run from within the 'tutorial' folder
      or the import statements need to be adjusted
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# add parent dir to path
p2parent = os.path.dirname( os.getcwd())
sys.path.append( p2parent)
from fractalAnalysis import synFractal, fractalDim

#=================================0==================
#                set parameters
#====================================================
geometry = 'koch' # toggle between koch and sierpinski
dPar = {  #----create synthetic fractal--------
          'length'    : 105, # length of fractal
          # --- C(r) and D2----------
          'x_min'     : 0.1, #smallest distance for point density calculation and C(r)
          #'Nmin'  : 3, # determine lower cut-off for fractal fit
          #-------plotting-----------------------
          'showFractal' : False,
          'plotFormat'  : 'png',
        }
if geometry == 'koch':
    dPar['a_nPoints'] = np.arange(3, 8, 1)
    dPar['D_true']    = np.log(4)/np.log(3)
elif geometry == 'sierpinski':
    dPar['a_nPoints'] = np.arange( 50, 10000, 1000)
    dPar['D_true'] = np.log(3)/np.log(2),
#=================================1=========================================================
#                             create fractal point clouds
#===========================================================================================
aN_tot = np.zeros( dPar['a_nPoints'].shape[0])
aD2    = np.zeros( dPar['a_nPoints'].shape[0])
iP = 0
aCol = plt.cm.get_cmap( plt.cm.RdYlGn)
for i, iN in enumerate( dPar['a_nPoints']):
    #vCantor = fractal.cantor(dPar['nPoints'])*dPar['length']
    if geometry == 'sierpinski':
        a_X, a_Y = synFractal.sierpinski( iN, length = dPar['length'],
                                    interactivePlot = False,
                                    plotTriangle = True)
    elif geometry == 'koch':
        a_X, a_Y = synFractal.kochCoastline( iN, length = dPar['length'])
    if dPar['showFractal']:
        plt.figure(10)
        plt.title( f"N={iN}")
        plt.plot( a_X, a_Y, 'ko', ms = 1)
        plt.axis( 'equal')
        #plt.pause( .5)
        plt.show()

    #=================================2=========================================================
    #                    Correlation Integral
    #===========================================================================================
    d_C_int = fractalDim.C_r_Int( a_X, a_Y,
                                 x_min = dPar['x_min'],
                                 x_max = dPar['length']*3,
                                logBinning = True)
    N_max = ( a_X.shape[0]*(a_X.shape[0]-1))*.5
    # determine max. cut-off
    dRmax = { 'rmax' : d_C_int['aL_ave'][d_C_int['aN_sum'] > N_max][0]}
    print( 'Nmax:', N_max, dRmax)
    sel_pl   = d_C_int['aL_ave'] < dRmax['rmax']
    #=================================3=========================================================
    #                        fractal dimension
    #===========================================================================================
    slope, interc, r, p, __ = scipy.stats.linregress( np.log10( d_C_int['aL_ave'][sel_pl]),
                                                      np.log10( d_C_int['aCorr'][sel_pl]))
    interc = 10**interc
    print( 'D2' , slope, 'rmax', dRmax['rmax'], f"true D ({geometry})", dPar['D_true'])
    aC_hat = interc*d_C_int['aL_ave']**slope
    #aC_hat = slope*np.log10( d_C_int['aL_ave']) + interc
    # store D2 for current number of points
    aN_tot[i] = a_X.shape[0]/dPar['length']
    aD2[i]    = slope
    #
    color  = aCol( i/len( dPar['a_nPoints']))
    plt.figure(1)
    ax = plt.subplot( 111)
    ax.loglog( d_C_int['aL_ave'], d_C_int['aCorr'], 'o',
                ms = 5, mec = color, mfc = 'w', label = f"N={iN}")
    ax.loglog( d_C_int['aL_ave'], aC_hat, '--', color = color)
# mark upper bound
ax.axvline( dRmax['rmax'])
#------labels and legends
ax.legend( loc = 'lower right')
ax.set_xlabel( 'Distance ()')
ax.set_ylabel( 'Correlation Integral')
plt.savefig( f"{p2parent}/plots/CorrInt_{geometry}.{dPar['plotFormat']}")

#
plt.figure(2)
ax2 = plt.subplot(111)
ax2.plot( aN_tot, aD2, 'ko', ms = 6, mfc = 'w', mew = 1.5, label = '$D_\mathrm{obs}$')
ax2.plot( ax2.get_xlim(), [dPar['D_true'], dPar['D_true']], 'r--', label = '$D_\mathrm{True}$')
ax2.legend( loc = 'lower right', frameon = False, title = None)
ax2.set_xlabel( '$\#$Ev./mm')
ax2.set_ylabel( 'Fractal Dimension')
ax2.set_ylim( 1, 2)
plt.savefig( f"{p2parent}/plots/D2_{geometry}.{dPar['plotFormat']}")
plt.show()







    
    