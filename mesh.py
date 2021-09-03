# Implement class mesh
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy import spatial
from tqdm import tqdm
from point2D import Point2D, compute_distance
from primitive_objects import Rectangle, Circle


class Mesh():
    """ Mesh class
    """

    def __init__(self):
        self.cells = []  # represents a collection of points of the mesh
        self.an = 0.0  # resolution/pi of the mesh (equivalent circle)
        # key is the pos in the vector of primitiveObj, value is idObj
        self.idBackground = {0: -1}
    #--------------------------------------------


    def __len__(self):
        return len(self.cells)
    #--------------------------------------------

    def __buildKdTree(self):
        # Build the kd-tree for all mid-points of the mesh
        numCells = len(self.cells)
        cellCoords = np.zeros((numCells, 2))
        print("Building the kd-tree for the mesh. . .")
        for i in tqdm(range(numCells)):
            cellCoords[i, :] = np.array([self.cells[i].x, self.cells[i].y])
            # build the tree
        self.tree = spatial.KDTree(cellCoords)
    #--------------------------------------------

    # Query if a point is within a cell
    def pointInCell(self, robs):
        """ Given an observation point, returns the closest cell and the closest distance
            to that tet by using the kd-tree

            Output:
                distance to the closest cell
                closest cell elem
        """
        point = np.zeros((2))
        point[0] = robs.x
        point[1] = robs.y
        d_cell = self.tree.query(point)
        return d_cell[0], d_cell[1]
    #--------------------------------------------

    def plotMesh(self):
        # Extract x and y from cells
        if not self.cells:
            raise RuntimeError("Empty cells")
        x = []
        y = []
        for i in range(len(self.cells)):
            x.append(self.cells[i].x)
            y.append(self.cells[i].y)
        plt.plot(x, y, 'k+')
        plt.show()
    #--------------------------------------------


    def __createMeshFromBackgroundRectangle(self, rectangle, resolution):
        if resolution == 0.0:
            raise ValueError("Resolution cannot be zero")
        self.an = resolution
        # We assume that rectangle is envoloped within an equivalent square
        xmin = rectangle.pMin.x
        xmax = rectangle.pMax.x
        ymin = rectangle.pMin.y
        ymax = rectangle.pMax.y
        # Compute the # of cells in x and y
        self.nx = int(round((xmax - xmin + resolution) / resolution))
        self.ny = int(round((ymax - ymin + resolution) / resolution))
        # Create the grid points
        x = np.linspace(xmin, xmax, self.nx)
        y = np.linspace(ymin, ymax, self.ny)
        # Loop and set the cells
        for i in range(self.nx):
            for j in range(self.ny):
                aPt = Point2D(x[i], y[j])
                self.cells.append(aPt)
    #--------------------------------------------

    def generateMesh(self, vecObj, resolution):
        """ This function generate the mesh given a list of primitive object
        """
        self.an = resolution / pi
        if not vecObj:
            raise RuntimeError("List of primitive objects is empty")
        print('Found {} primitive objects'.format(len(vecObj)))
        # Extract the id for finding the minimu later
        minId = {}
        for i in range(len(vecObj)):
            try:
                minId[i].append(vecObj[i].getId())
            except KeyError:
                minId[i] = [vecObj[i].getId()]
        keyMin = min(minId, key=minId.get)
        self.idBackground = {keyMin: minId[keyMin]}
        # Generate mesh
        if type(vecObj[keyMin]) is Rectangle:
            self.__createMeshFromBackgroundRectangle(
                vecObj[keyMin], resolution)
        else:
            raise RuntimeError("Wrong background object")
        # Build the tree
        self.__buildKdTree()
    #--------------------------------------------

    def __str__(self):
        """ Printing the mesh
        """
        sMesh = 'Mesh:'
        sMesh += "\n"
        sMesh += 'Nx: {}'.format(self.nx)
        sMesh += "\n"
        sMesh += 'Ny: {}'.format(self.ny)
        sMesh += "\n"
        for i in range(len(self.cells)):
            sMesh += '{}: {}'.format(i, self.cells[i])
            sMesh += "\n"
        return sMesh
    #--------------------------------------------

# Quickly test
if __name__ == '__main__':

    rec1 = Rectangle(idObj=2)

    pMin = Point2D(0, 0)
    pMax = Point2D(2, 2)
    rec2 = Rectangle(idObj=1, pMin=pMin, pMax=pMax)
    circ1 = Circle(idObj=4)
    circ2 = Circle(idObj=3, radius=0.5, center=Point2D(0.0, 0.0))
    # append a couple of mix objs
    print("\n")
    vPrObj = []
    vPrObj.append(rec1)
    vPrObj.append(rec2)
    vPrObj.append(circ1)
    vPrObj.append(circ2)
    # Create mesj
    msh = Mesh()
    resolution = 0.25
    msh.generateMesh(vPrObj, resolution)
    msh.plotMesh()
    print(msh.pointInCell(Point2D()))
    print('No. of cells in the mesh: {}'.format(len(msh)))
    print(msh)