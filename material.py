import numpy as np
from math import pi
from object import Object
from mesh import Mesh
from em_utils import eps0


class Material(Object):
    """ A class that implements a dielectric material
    """

    def __init__(self, idObj=-1, epsr=1.0, sigma=0.0):
        Object.__init__(self, idObj)
        self.epsr = epsr
        self.sigma = sigma
        self.epsrc = 0.0 + 1j*0.0

    def setPermittivity(self, freq):
        """ Set the permittivity
        """
        if freq <= 0.0:
            raise ValueError("Frequency must be greater than zero")
        omega = 2.0*pi*freq
        self.epsrc = np.complex(self.epsr, -self.sigma/(omega*eps0))

    def __str__(self):
        """print info about the material
        """
        s0 = 'Material with idObj: {}'.format(self.getId())
        s2 = '(epsr, sigma): {}; {}'.format(self.epsr, self.sigma)
        s3 = 'epsr complex: {}'.format(self.epsrc)
        return s0 + "\n" + s2 + "\n" + s3


def set_dielectrics(vecMat, vecFreq, matToFreqMap):
    """ Modify the materials in order to set the complex permittivity
    """
    numMat = len(vecMat)
    numFreq = len(vecFreq)
    for idMat, idFreq in matToFreqMap.items():
        # Find this freq in the vec of freq
        for f in range(numFreq):
            if idFreq == vecFreq[f].getId():
                freq = vecFreq[f].freq
        # Set this mat
        for m in range(numMat):
            if idMat == vecMat[m].getId():
                vecMat[m].setPermittivity(freq)


def set_dielectric_properties_at_cells(msh, vecObj, vecMat, objToMatMap):
    """ Set the dielectric properties for each cell of the domain
    """
    numCells = len(msh.cells)
    numObj = len(vecObj)
    numMat = len(vecMat)
    dielectric = np.zeros((numCells), dtype=np.complex)
    # Start loop over the cell
    for i in range(numCells):
        foundCell = False
        for j in range(numObj):
            if foundCell == True:
                break
            if vecObj[j].isInternal(msh.cells[i]) == True:
                idMat = objToMatMap[vecObj[j].getId()]
                for m in range(numMat):
                    if idMat == vecMat[m].getId():
                        dielectric[i] = vecMat[m].epsrc
                        foundCell = True
                        break
    return dielectric


if __name__ == '__main__':
    mat = Material()
    print(mat)
    mat.setPermittivity(1e9)
    print(mat)
