import copy
import numpy as np
from tqdm import tqdm
from scipy import special 
from scipy import fft
from scipy.sparse.linalg import LinearOperator
from math import pi
from point2D import compute_distance
from em_utils import get_dielectric_wavenumber
from material import compute_contrast

class InhomogeneousCylinderKernel(LinearOperator):
    """ A class for TM Wave Scattering from Inhomogeneous Dielectric Cylinders: Volume EFIE
        The Integral Kernel is: A = I + Z*Chi
            where Chi=diag(chi), and I is diag(1)
    """

    def __init__(self, msh, freq, epsrcBackground, epsrc, assembleFullZ=False):
        """ Init the kernel
        """
        self.dofs = len(msh)
        self.nx = msh.nx
        self.ny = msh.ny

        # Need to extend this from LinearOperator
        super().__init__(np.complex, (self.dofs, self.dofs))
        
        self.kb = get_dielectric_wavenumber(freq, epsrcBackground)  # k of background
        self.chi = compute_contrast(epsrcBackground, epsrc) # contrast vector

        # Compute and store the FFT of the circulant Z-matrix
        Zext = self.__assembleZextended(msh)
        self.__Zpf = fft.fft2(Zext) # keep the FFT of the extedend matrix

        # In case requested, assemble the full Z matrix
        if assembleFullZ:
            self.Zmat = self.__assembleZmat(msh)

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
            return (1j*pi*kb*an*0.5)*special.jv(1,kb*an)*special.hankel2(0,kb*dist)
    #--------------------------------------------

    def assembleZmatRow(self, irow, msh):
        """ Assemble one row of the impedance matrix
        """

        an = msh.an # get the resolution
        Zrow = np.zeros((self.getDoFs()), dtype=np.complex)
        for jcol in range(self.getDoFs()):
            Zrow[jcol] = self.__computeZmn(self.kb, an, irow, jcol, msh.cells[irow], msh.cells[jcol])

        return Zrow
    #--------------------------------------------

    def __assembleZmat(self, msh):
        """ Assemble the impedance matrix (full)
        """
        Zmat = np.zeros((self.getDoFs(), self.getDoFs()), dtype=np.complex)
        print('Assembling the impedance matrix. . .')
        for irow in tqdm(range(self.getDoFs())):
            Zmat[irow, :] = self.assembleZmatRow(irow, msh)
        # Retunr the matrix
        return Zmat
    #--------------------------------------------

    def __assembleZextended(self, msh):
        """ Assemble the Z matrix in a circulant manner
        """
        Zrow = self.assembleZmatRow(0, msh)  # first compute a row
        return self.__extendedZmat(Zrow, self.nx, self.ny)
    #--------------------------------------------

    def getDoFs(self):
        return self.dofs
    #--------------------------------------------

    def __extendedIndex(self, m, Nm):
        # Calculate the p'index in the circulant matrix
        if m >= 0 and m <= Nm-1:
            return m
        else:
            return 2*Nm - (m+1)
    #--------------------------------------------

    def __extendedZmat(self, Zrow_, Nx, Ny):
        Px = 2*Nx-1
        Py = 2*Ny-1
        Zp = np.zeros((Px, Py), dtype=np.complex)
        Zrow = Zrow_.reshape(Nx, Ny)
        for p in range(Px):
            # Calculate the p'index in the circulant matrix
            pp = self.__extendedIndex(p, Nx)
            for q in range(Py):
                # Calculate the q' index in the circulant matrix
                qq = self.__extendedIndex(q, Ny)

                # Populate the circulant Z
                Zp[p, q] = Zrow[pp, qq]
        # Return the circ matrix
        return Zp
    #--------------------------------------------

    def __extendedVector(self, v_, Nx, Ny):
        """ Create an extended vector in a circulant fashion
        """
        Px = 2*Nx-1
        Py = 2*Ny-1
        vp = np.zeros((Px, Py), dtype=np.complex)
        v = v_.reshape(Nx, Ny)
        for p in range(Px):
            for q in range(Py):
                if (p >= 0 and p <= Nx-1) and (q >= 0 and q <= Ny-1):
                    vp[p, q] = v[p, q]
        # return the vector
        return vp
    #--------------------------------------------

    def _matvec(self, x):
        """ This function overloads the matrix-vector product.
            It implements a FAST A*x operations using the FFT.
            A = I + Z*Chi, then A*x = x + Z*p, where p = chi .* x
        """
        # Vector FFT operations
        if x.shape != self.chi.shape:
            raise RuntimeError("_matvec DIO PORCO")

        p = np.multiply(self.chi, x)
        vp = self.__extendedVector(p, self.nx, self.ny) # circulant vector
        # Forward FFT
        vpf = fft.fft2(vp)

        # Hadamard prod of the FFTs
        zvf = np.multiply(self.__Zpf,vpf)
        # Inverse FFT
        izvf = fft.ifft2(zvf)
        # Get the vector back from the circulant matrix
        xzf = izvf[0:self.nx, 0:self.ny]
        zl = xzf.reshape(self.nx*self.ny)
        return x + zl
    #--------------------------------------------

    def full_matvec(self, x):
        """ Compute the mat-vec product using the full matrix
            Use this function only for testing
        """
        if hasattr(self, 'Zmat'):
            p = np.multiply(self.chi, x)
            Zp = self.Zmat.dot(p)
            return x + Zp
        else:
            raise Warning("Mat-Vec with Full Zmat: Missing Zmat")
    #--------------------------------------------
#-----------------------------------------------------------------------------------------