import copy
import numpy as np
from tqdm import tqdm
from scipy import special
from scipy import fft 
from scipy.sparse.linalg import LinearOperator
from math import pi
from point2D import compute_distance
from em_utils import get_impedance, get_k0

class InhomogeneousCylinderKernel(LinearOperator):
    """ A class for TM-Wave Scattering from Inhomogeneous Dielectric Cylinders: Volume EFIE 
    """

    def __init__(self, msh, freq, dielProp, dtype=np.complex):
        """ Init the kernel
        """

        self.dofs = len(msh)
       
        # Need to extend this from LinearOperator
        self.explicit = False
        self.shape = (self.dofs, self.dofs)
        self.dtype = dtype

        self.eta = get_impedance(epsrc=1.0)  # background is air / free space
        self.kb = get_k0(freq) # k of background
        self.freq = freq
        self.prop = copy.copy(dielProp)  # build a swallow copy
    #--------------------------------------------
    
    def __computeZmn(self, kb, an, m, n, cm, cn):
        """ Compute an entry for the Z matrix.
            It can be self or off-diagonal entries depending on (m,n).
            This function is private.
        """
        if m==n: #self term
            return (1j*0.5)*(pi*an*kb*special.hankel2(1, kb*an) -2*1j)
        else: # off-diagonal
            dist = compute_distance(cm, cn)
            return (1j*pi*kb*an*0.5)*special.j1(kb*an)*special.hankel2(0,kb*dist)
    #--------------------------------------------


    def assembleZmatRow(self, irow, msh):
        """ Assemble one row of the impedance matrix
        """
        # an = 
        # kb = 

        Zrow = np.zeros((self.getDoFs()), dtype=np.complex)
        for jcol in range(self.getDoFs()):
            Zrow[jcol] = self.__computeZmn(kb, an, irow, jcol, msh.cells[irow], msh.cells[jcol])

        return Zrow
    #--------------------------------------------

    def assembleZmat(self, msh):
        """ Assemble the impedance matrix
        """
        Zmat = np.zeros((self.getDoFs(), self.getDoFs()), dtype=np.complex)
        print('Assembling the impedance matrix. . .')
        for irow in tqdm(range(self.getDoFs())):
            Zmat[irow, :] = self.assembleZmatRow(irow, msh)
        # Retunr the matrix
        return Zmat
    #--------------------------------------------

    
    def getDoFs(self):
        return self.dofs
    #--------------------------------------------

