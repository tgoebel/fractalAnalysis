# !bin/python3.7
import matplotlib.pyplot as plt
from numpy import array, append

def kochenize(a,b):
    HFACTOR = (3**0.5)/6
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    mid = ( (a[0]+b[0])/2, (a[1]+b[1])/2 )
    p1 = ( a[0]+dx/3, a[1]+dy/3 )
    p3 = ( b[0]-dx/3, b[1]-dy/3 )
    p2 = ( mid[0]-dy*HFACTOR, mid[1]+dx*HFACTOR )
    return p1, p2, p3

def koch( steps, width, verbose = False, **kwargs):
    """ half snowflake: fractal dimension: log4/log3 """
    arraysize = 4**steps + 1
    points = [(0.0,0.0)]*arraysize
    points[0] = (-width/2., 0.0)
    points[-1] = (width/2., 0.0)
    stepwidth = arraysize - 1
    vX, vY = array([]), array([])
    if 'interactivePlot' in kwargs.keys() and kwargs['interactivePlot'] == True:
        plt.figure( 1, figsize = (12,5))
        plt.clf()
        ax = plt.axes( [.1,.1,.85,.85]) #[0, 10, 0, 1])
        plt.ion()


    for n in range( steps):
        segment = int( (arraysize-1)/stepwidth)
        if verbose:
            print('steps', steps, 'segment', segment)

        for s in range( segment):
            st = int(s*stepwidth)
            a = (points[st][0], points[st][1])
            b = (points[st+stepwidth][0], points[st+stepwidth][1])
            index1 = int( st + (stepwidth)/4)
            index2 = int( st + (stepwidth)/2)
            index3 = int( st + ((stepwidth)/4)*3)
            result = kochenize(a,b)
            vX = append( vX, points[st][0])
            vY = append( vY, points[st][1])
            points[index1], points[index2], points[index3] = result    
        if 'interactivePlot' in kwargs.keys() and kwargs['interactivePlot'] == True:
            ax.set_title(  'N = %i'%(vX.shape[0]))
            ax.plot( vX, vY, 'ko', ms = 1)
            #ax.text( B[0], C[1],  'N = %i'%(n))
            plt.pause(0.05)
        stepwidth = int( stepwidth/4)
    return vX, vY

if __name__ == '__main__':
    TOTALWIDTH = 10
    vX, vY = koch( 5, TOTALWIDTH )
    plt.plot( vX, vY, 'ko', ms=1, )
    plt.axis('equal')
    plt.show()




