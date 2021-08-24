from tqdm import tqdm
from scipy import special
from scipy import fft
import numpy as np
from numpy import linalg as LA
from math import pi, sqrt
import matplotlib.pyplot as plt


class Cell:
    """ A cell is a point with x and y coordinates
    """
    def __init__(self, rx, ry):
        self.rx = rx
        self.ry = ry

    # Print on the class
    def __str__ (self):
        return "cell(x=" + str(self.rx) + ", y=" + str(self.ry) + ")" + "\n"
#--------------------------------------------


def distCells(cm, cn):
    """ Compute the distance between two cells
    """
    return sqrt((cm.rx-cn.rx)**2 + (cm.ry-cn.ry)**2)
#--------------------------------------------


def testMesh(xmin, xmax, ymin, ymax, Nx, Ny):
    xCoord = np.linspace(xmin, xmax, Nx)
    yCoord = np.linspace(ymin, ymax, Ny)
    dx = xCoord[1] - xCoord[0]
    dy = yCoord[1] - yCoord[0]
    Sm = dx*dy
    resolution = sqrt(Sm)/pi

    # Set an array of cells
    cells = np.zeros((Nx*Ny), dtype=Cell)
    cnt = 0
    for i in range(Nx):
        for j in range(Ny):
            cells[cnt] = Cell(xCoord[i], yCoord[j])
            cnt = cnt + 1

    return resolution, cells
#--------------------------------------------


def computeZmn(kb, an, m, n, cm, cn):
    """ Compute an entry for the Z matrix.
        It can be self or off-diagonal entries depending on (m,n)
    """
    if m==n: #self term
        return (1j*0.5)*(pi*an*kb*special.hankel2(1, kb*an) -2*1j)
    else: # off-diagonal
        dist = distCells(cm, cn)
        return (1j*pi*kb*an*0.5)*special.j1(kb*an)*special.hankel2(0,kb*dist)

def computeRowZ(kb, an, irow, cells):
    """ Compute a row of the Z-matrix
    """
    N = len(cells)
    zrow = np.zeros((N), dtype=np.complex)
    for jcol in range(N):
        zrow[jcol] = computeZmn(kb, an, irow, jcol, cells[irow], cells[jcol])
    return zrow
#--------------------------------------------

def computeFullZ(kb, an, cells):
    N = len(cells)
    Zmat = np.zeros((N,N), dtype=np.complex)
    # Loop over the row and assemble the matrix
    for irow in tqdm(range(N)):
        Zmat[irow, :] = computeRowZ(kb, an, irow, cells)
    return Zmat
#--------------------------------------------


def extendedIndex(m, Nm):
    # Calculate the p'index in the circulant matrix
    if m>=0 and m<=Nm-1:
        return m
    else:
        return 2*Nm - (m+1)

#--------------------------------------------


def extendedZmat(Zrow_, Nx, Ny):
    Px = 2*Nx-1
    Py = 2*Ny-1
    Zp = np.zeros((Px, Py), dtype=np.complex)
    Zrow = Zrow_.reshape(Nx, Ny)
    for p in range(Px):
        # Calculate the p'index in the circulant matrix
        pp = extendedIndex(p, Nx)
        for q in range(Py):
            # Calculate the q' index in the circulant matrix
            qq = extendedIndex(q, Ny)

            # Populate the circulant Z
            Zp[p,q] = Zrow[pp,qq]
    # Return the circ matrix
    return Zp
#--------------------------------------------

def extendedVector(v_, Nx, Ny):
    """ Create an extended vector in a circulant fashion
    """
    Px = 2*Nx-1
    Py = 2*Ny-1
    vp = np.zeros((Px, Py), dtype=np.complex)
    v = v_.reshape(Nx, Ny)
    for p in range(Px):
        for q in range(Py):
            if (p>=0 and p <=Nx-1) and (q>=0 and q <=Ny-1):
                vp[p,q] = v[p,q]

    # return the vector
    return vp
#--------------------------------------------


# Test
if __name__ == '__main__':
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    N = 20
    an, cells = testMesh(xmin, xmax, ymin, ymax, N, N)
    # for i in range(N*N):
    #     print(cells[i])
    print("resolution: ", distCells(cells[3], cells[0]))

    # Compute one row of the matrix
    freq = 1e9
    lam = 299792458.0 / freq
    kb = (2*pi)/lam
    z0 = computeRowZ(kb, an, 0, cells)
    # print(z0)
    z1 = computeRowZ(kb, an, 1, cells)
    # print(z1)
    z3 = computeRowZ(kb, an, 3, cells)
    # print(z3)
    print(LA.norm(z0[3]-z3[0]))
    print(LA.norm(z0[1]-z1[0]))
    print(LA.norm(z1[3]-z3[1]))

    # Z full and vector
    v = np.ones((N*N), dtype=np.complex)
    Zmat = computeFullZ(kb, an, cells)
    # plt.matshow(np.abs(Zmat))
    x = Zmat.dot(v)

    # Circulant matrix and vector
    Zp = extendedZmat(z0, N, N) #use the first row
    vp = extendedVector(v, N, N)

    Zpf = fft.fft2(Zp)
    vpf = fft.fft2(vp)
    zvf = np.multiply(Zpf, vpf)
    izvf = fft.ifft2(zvf)
    xzf = izvf[0:N, 0:N]
    xfast = xzf.reshape(N*N)
    # Norm fast and full
    print("Norm fast-full: ", LA.norm(xfast-x))
