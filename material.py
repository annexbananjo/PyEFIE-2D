import numpy as np
import matplotlib.pyplot as plt
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
#--------------------------------------------


def __set_dielectrics(vecMat, vecFreq, matToFreqMap):
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
#--------------------------------------------


def set_material_id_at_cells(msh, vecObj, vecMat, objToMatMap):
    """ At each cell of the domain, set the id of the material associated to 
        that cell
    """
    numCells = len(msh.cells)
    numObj = len(vecObj)
    numMat = len(vecMat)
    dielectric = np.zeros((numCells), dtype=np.complex)

    # create an array with the ID of the material
    idBackground = msh.idBackground[0][0]
    vecMatID = idBackground * np.ones((numCells), dtype=np.int8)

    # Start loop over the cell
    for i in range(numCells):
        foundCell = False
        for j in range(numObj):
            if foundCell == True:
                break
            # skip the background obj
            if vecObj[j].idObj != idBackground:
                if vecObj[j].isInternal(msh.cells[i]) == True:
                    idMat = objToMatMap[vecObj[j].getId()]
                    vecMatID[i] = idMat
                    foundCell = True # found only one time
    return vecMatID
#--------------------------------------------


def set_dielectric_properties_at_cells(vecFreq, matToFreqMap, vecMat, vecMatID):
    """ Set the dielectric properties at each cell
    """

    # First set the dielectric properties for the material database
    __set_dielectrics(vecMat, vecFreq, matToFreqMap)

    numCells = len(vecMatID)
    numMat = len(vecMat)
    epsrc = np.zeros((numCells), dtype=np.complex)
    for i in range(numCells):
        idMat = vecMatID[i]
        foundDiel = False
        for m in range(numMat):
            if idMat == vecMat[m].getId():
                epsrc[i] = vecMat[m].epsrc
                foundDiel = True
        if foundDiel==False:
            raise ValueError("Cannot find dielectric in the material database")
    return epsrc
#--------------------------------------------
    

def plot_cells_with_ID(msh, vecMatID):
        """ Given the mesh and a vector with size numCells, 
            print the mesh with the material ID associated to each cell
        """
        numCells = len(msh.cells)
        if len(vecMatID) != numCells:
            raise ValueError("vecMatID is not coherent with the mesh")
        x = []
        y = []
        txt = []
        for i in range(len(msh.cells)):
            x.append(msh.cells[i].x)
            y.append(msh.cells[i].y)
            txt.append(vecMatID[i])
        plt.plot(x, y, '.')

        # zip joins x and y coordinates in pairs
        cnt = 0
        for xs,ys in zip(x,y):
            plt.annotate(txt[cnt], # this is the text
                        (xs,ys), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(1,1), # distance from text to points (x,y)
                        ha='center') # horizontal alignment can be left, right or center
            cnt+=1
        plt.show()
#--------------------------------------------

def get_epsrc_background(msh, vecMat):
    """ Gets the epsrc of the background
    """
    numMat = len(vecMat)
    for m in range(numMat):
        if msh.idBackground[0][0] == vecMat[m].getId():
            epsrc = vecMat[m].epsrc
    return epsrc
#--------------------------------------------


def compute_contrast(epsrcBackground, epsrc):
    """ Compute the contrast chi = (epsrc(r) - epsrcBackgroung) / epsrcBackground
    """
    N = len(epsrc)
    chi = np.zeros((N), dtype=np.complex)
    for i in range(N):
        chi[i] = (epsrc[i] - epsrcBackground)/ epsrcBackground
    return chi
#-----------------------------------------------------------------------------------------


if __name__ == '__main__':
    mat = Material()
    print(mat)
    mat.setPermittivity(freq=1e9)
    print(mat)
