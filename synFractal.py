#!/bin/python3.7
"""
    - methods to compute synthetic fractals with known dimensions
      1) cantor
      2) koch coastline
      3) sierpinski triangles

"""
import numpy as np
import matplotlib.pyplot as plt
# ===============================================================================
#                     methods to create synthetic fractals
# ===============================================================================
def kochCoastline( nSteps, **kwargs):
    """ koch coastline: fractal dimension: log4/log3 = 1.26185
    :Input - nIter  =  int,  --> iterations not actual number of points < 10
            kwargs:
            'length'   : float , default = 10
    :return  vX, vY - 2D coordinates of points creating coastline
    """
    length = 100
    if 'length' in kwargs.keys() and kwargs['length'] is not None:
        length = kwargs['length']
    from fractalAnalysis import kochCoastline as koch
    vX, vY = koch.koch( nSteps, length, **kwargs)  # vX+abs(vX.min()), vY
    return vX, vY

def sierpinski(nPoints, **kwargs):
    """ sierpinski triangles: fractal dimension: log3/log2 = 1.5849
    :Input - nPoints : int #(number of iterations)
                kwargs:
                     'length'         : float
                     'plotTriangle'   : True
                     interactivePlot' : True  # animation
               'startingPoint'        : int,

    """
    startingPoint = 0
    if 'startingPoint' in kwargs.keys() and kwargs['startingPoint'] is not None:
        startingPoint = kwargs['startingPoint']
    length = 10
    if 'length' in kwargs.keys() and kwargs['length'] is not None:
        length = kwargs['length']
    h = (.5 * np.sqrt(3) * length)
    # matrix for dot product which provides coordinates of new points
    halfVec = np.array([[.5, startingPoint],
                        [startingPoint,.5]]);  # go half distance means dot product with vector and standard unit vector time .5
    A = np.array([startingPoint, startingPoint])
    B = np.array([length, startingPoint])
    C = np.array([length / 2, h])
    mCorners = np.vstack((A, B, C))

    # starting point
    x = mCorners[startingPoint - 1]
    # Note the comma after line. This is placed here because plot returns a list of lines that are drawn.
    if 'interactivePlot' in kwargs.keys() and kwargs['interactivePlot'] == True:
        plt.clf()
        ax = plt.axes([.1, .1, .85, .85])  # [0, 10, 0, 1])
        plt.ion()
        if kwargs['plotTriangle']:
            ax.plot([A[0], B[0]], [A[1], B[1]], 'ko-')
            ax.plot([B[0], C[0]], [B[1], C[1]], 'ko-')
            ax.plot([C[0], A[0]], [C[1], A[1]], 'ko-')
            #line, = ax.plot(x[0], x[1], 'ko', markersize=2)

    vX = np.array([])
    vY = np.array([])
    # start iteration
    n = 1
    while n < nPoints:
        corner = np.random.randint(1, 4)
        # new point is HALF of the (randomly chosen corner - coordinates of old points) + coordinates of old point
        x = np.dot(halfVec, (mCorners[corner - 1] - x)) + x
        vX = np.append(vX, x[0])
        vY = np.append(vY, x[1])
        if kwargs['interactivePlot'] == True and n % 100 == 0:
            # plt.scatter(i, y)
            ax.set_title('N = %i' % (n))
            ax.plot(vX, vY, 'ko', ms=1)
            # ax.text( B[0], C[1],  'N = %i'%(n))
            plt.pause(0.05)
        # print( '#n iterations',n+1
        n += 1
    return vX, vY


def cantor(n):
    """ create 1D cantor set:
    D = ln2/ln3 = 0.6309
     """
    return np.array([0.] + cant(0., 1., n) + [1.])


def cant(x, y, n):
    if n == 0:
        return []

    new_pts = [2. * x / 3. + y / 3., x / 3. + 2. * y / 3.]
    return cant(x, new_pts[0], n - 1) + new_pts + cant(new_pts[1], y, n - 1)
