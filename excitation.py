import numpy as np
from math import cos, sin
from cmath import exp

class Planewave():
    """ Class that handles a simple planewave
    """
    def __init__(self, numPlws, startPhi, endPhi):
        """ Init the object planewave
        """
        self.numPlws = numPlws
        self.startPhi = np.deg2rad(startPhi)
        self.endPhi = np.deg2rad(endPhi)
        self.phiInc = np.linspace(startPhi, endPhi, numPlws)
    #--------------------------------------------

    def computeEinc(self, kb, msh):
        """ Compute the incident electric field at each point of the mesh
        """
        if self.numPlws <=0:
            raise RuntimeError("No. of excitations (planewave) must be > 0")

        numCells = len(msh)
        Einc = np.zeros((numCells,self.numPlws), dtype=np.complex)

        # start looping for each cell and plw
        for i in range(self.numPlws):
            cosPhi = np.cos(self.phiInc[i])
            sinPhi = np.sin(self.phiInc[i])
            for j in range(numCells):
                xc = msh.cells[j].x * cosPhi
                yc = msh.cells[j].y * sinPhi
                e0 = exp(-1j*kb*(xc + yc))
                Einc[j,i] =e0
        return Einc
    #--------------------------------------------



