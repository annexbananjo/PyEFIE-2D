import copy
import numpy as np
from tqdm import tqdm
from scipy.special import hankel2, jv
from math import pi
from point2D import compute_distance
from em_utils import get_impedance, get_k0


class IEKernel():
    """ Abstract class
    """

    def __init__(self, msh):
        """ Init the kernel with mesh and the dielectric or contrast prop
        """
        self.dofs = len(msh)

    def getDoFs(self):
        return self.dofs

    def assembleZmatRow(self, irow, msh):
        raise NotImplementedError()

    def assembleZmat(self, msh):
        raise NotImplementedError()


class InhomogeneousCylinderKernel(IEKernel):
    """ A class for TM-Wave Scattering from Inhomogeneous Dielectric Cylinders: Volume EFIE 
    """

    def __init__(self, msh, freq, dielProp):
        """ Ini the kernel
        """
        IEKernel.__init__(self, msh)

        self.eta = get_impedance(1.0)  # background is air / free space
        self.k = get_k0(freq)
        self.coef_mn = self.eta * pi * msh.an * 0.5 * jv(1, self.k * msh.an)
        self.coef_mm = self.eta * pi * msh.an * \
            0.5 * hankel2(1, self.k * msh.an)
        self.freq = freq
        self.prop = copy.copy(dielProp)  # build a swallow copy

    def assembleZmatRow(self, irow, msh):
        """ Assemble one row of the impedance matrix"""
        Zrow = np.zeros((self.getDoFs()), dtype=np.complex)
        thisCell = msh.cells[irow]
        # Start loop over the cells
        for i in range(self.getDoFs()):
            if i == irow:
                # self term
                Zrow[i] = self.coef_mm - 1j*self.eta * \
                    self.prop[i] / (self.k * (self.prop[i] - 1.0))
            else:
                # off-diag
                Zrow[i] = self.coef_mn * \
                    hankel2(0, compute_distance(thisCell, msh.cells[i]))
        return Zrow

    def assembleZmat(self, msh):
        """ Assemble the impedance matrix
        """
        Z = np.zeros((self.getDoFs(), self.getDoFs()), dtype=np.complex)
        print('Assembling the impedance matrix. . .')
        for i in tqdm(range(self.getDoFs())):
            Z[i, :] = self.assembleZmatRow(i, msh)
        # Retunr the matrix
        return Z
