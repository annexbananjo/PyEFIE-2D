import unittest
import numpy as np
from numpy import linalg as LA
from parse_json import parse_json
from mesh import Mesh
from object import Object
from point2D import Point2D
from material import Material, plot_cells_with_ID
from material import set_material_id_at_cells, set_dielectric_properties_at_cells
from material import get_epsrc_background
from frequency import Frequency
from em_utils import get_dielectric_wavenumber
from inhomogeneous_cylinder_kernel import InhomogeneousCylinderKernel

class TestPyEFIE2D(unittest.TestCase):
    # Test assembling rows
    def test_assembly_rows(self):
        # First parse the json file
        fileJSON = 'data/PetersonCylinder/PetersonCyl.json'
        resolution, vecObj, vecInc, vecFreq, vecMat, objToMatMap, matToFreqMap = parse_json(
            fileJSON)

        # Create the mesh
        msh = Mesh()
        msh.generateMesh(vecObj, resolution)

        # Set the dielectrics
        vecMatID = set_material_id_at_cells(msh, vecObj, vecMat, objToMatMap)
        epsrc = set_dielectric_properties_at_cells(vecFreq, matToFreqMap, vecMat, vecMatID)
        epsrcBackground = get_epsrc_background(msh, vecMat)

        # getting the frequency of the simulation
        freq = vecFreq[0].freq

        # Compute the Impedance matrix
        Zobj = InhomogeneousCylinderKernel(msh, freq, epsrcBackground, epsrc)

        # compare rows
        z0 = Zobj.assembleZmatRow(0, msh)
        z1 = Zobj.assembleZmatRow(1, msh)
        z3 = Zobj.assembleZmatRow(3, msh)

        tol = 1.0e-12
        self.assertLess(LA.norm(z0[3]-z3[0]), tol)
        self.assertLess(LA.norm(z0[1]-z1[0]), tol)
        self.assertLess(LA.norm(z1[3]-z3[1]), tol)
        #--------------------------------------------
    
    # Test assembling FULL vs. FAST
    def test_assembly_fast(self):
        # First parse the json file
        fileJSON = 'data/PetersonCylinder/PetersonCyl.json'
        resolution, vecObj, vecInc, vecFreq, vecMat, objToMatMap, matToFreqMap = parse_json(
            fileJSON)

        # Create the mesh
        msh = Mesh()
        msh.generateMesh(vecObj, resolution)

        # Set the dielectrics
        vecMatID = set_material_id_at_cells(msh, vecObj, vecMat, objToMatMap)
        epsrc = set_dielectric_properties_at_cells(vecFreq, matToFreqMap, vecMat, vecMatID)
        epsrcBackground = get_epsrc_background(msh, vecMat)

        # getting the frequency of the simulation
        freq = vecFreq[0].freq

        # Compute the Impedance matrix full
        Zobj = InhomogeneousCylinderKernel(msh, freq, epsrcBackground, epsrc, assembleFullZ=True)

        # Create a RHS
        x = np.ones((len(msh)),dtype=np.complex)
        y = Zobj.full_matvec(x)

        yy = Zobj.dot(x)

        tol = 1.0e-12
        self.assertLess(LA.norm(y-yy), tol)
        #--------------------------------------------

#--------------------------------------------

# Run the UT
if __name__ == '__main__':
    unittest.main()
       