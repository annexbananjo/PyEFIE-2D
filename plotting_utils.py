import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from point2D import Point2D


def plot_complex_quantity_on_2D_grid(msh, F, ngridx1=50, ngridx2=50):
    """
    """

    numCells = len(msh)
    x1 = []
    x2 = []
    for i in range(numCells):
        x1.append(msh.cells[i].x)
        x2.append(msh.cells[i].y)

    # Create grid values first.
    x1m = np.min(x1)
    x1p = np.max(x1)
    x2m = np.min(x2)
    x2p = np.max(x2)

    # Start plotting
    fig, (ax1,ax2) = plt.subplots(nrows=2)
    # -----------------------
    # Interpolation on a grid
    # -----------------------
    # A contour plot of irregularly spaced data coordinates
    # via interpolation on a grid.
    x1i = np.linspace(x1m, x1p, ngridx1)
    x2i = np.linspace(x2m, x2p, ngridx2)

    # Perform linear interpolation of the data (x,y)
    # on a grid defined by (xi,yi)
    triang = tri.Triangulation(x1, x2)
    interpolatorReal = tri.LinearTriInterpolator(triang, F.real)
    interpolatorIm = tri.LinearTriInterpolator(triang, F.imag)
    grid1, grid2 = np.meshgrid(x1i, x2i)
    fiReal = interpolatorReal(grid1, grid2)
    fiIm = interpolatorIm(grid1, grid2)

    # Note that scipy.interpolate provides means to interpolate data on a grid
    # as well. The following would be an alternative to the four lines above:
    #from scipy.interpolate import griddata
    #zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    #cmap="RdBu_r"
    ax1.contour(grid1, grid2, fiReal, levels=32, linewidths=0.5, colors='k')
    cntr1 = ax1.contourf(grid1, grid2, fiReal, levels=14, cmap="RdBu_r")

    ax2.contour(grid1, grid2, fiIm, levels=32, linewidths=0.5, colors='k')
    cntr2 = ax2.contourf(grid1, grid2, fiIm, levels=14, cmap="RdBu_r")

    fig.colorbar(cntr1, ax=ax1)
    ax1.set(xlim=(x1m, x1p), ylim=(x2m, x2p))
    strTitle = "Real(F)"

    ax1.set_title(strTitle)
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel("(m)")
    plt.ylabel("(m)")
    plt.show()

    fig.colorbar(cntr2, ax=ax2)
    ax2.set(xlim=(x1m, x1p), ylim=(x2m, x2p))
    strTitle = "Im(F)"

    ax2.set_title(strTitle)
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel("(m)")
    plt.ylabel("(m)")
    plt.show()
#-----------------------------------------------------------------------------------------

def plot_efield_in_points(xp, yp, msh, E):
    """ Plot the E-fields in the points (xp,yp)
    """
    numCells = len(msh)
    Np = len(xp)
    N_ = len(yp)
    if N_!=Np:
        raise RuntimeError("xp and yp must have same length")

    # Search in the list of points in the mesh
    Ei = []
    points = []
    for i in range(Np):
        pi = Point2D(xp[i], yp[i])
        points.append(pi)
        _,ci = msh.pointInCell(pi)
        Ei.append(E[ci])
    return points, Ei
