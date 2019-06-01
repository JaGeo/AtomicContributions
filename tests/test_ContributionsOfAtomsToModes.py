from ContributionsOfAtomsToModes import AtomicContributionToModes

import unittest
import math
import numpy as np

class AtomicContributionToModesTest(unittest.TestCase):
    def setUp(self):
        self.Contributions = AtomicContributionToModes(PoscarName='POSCAR', ForceConstants=False, ForceFileName='FORCE_SETS',
                                         supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],
                                         primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # self.Contributions2 = AtomicContributionToModes(PoscarName='POSCAR', ForceConstants=False, ForceFileName='FORCE_SETS',
        #                                  supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],nac='BORN',
        #                                  primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])


    def test_attributes(self):
        #test calculation of frequencies
        self.assertAlmostEqual(self.Contributions._frequencies[47], 3490.6434922723,places=1)
        #test calculation of eigenvectors
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 0]),0.00084433323436)
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 1]),0.00084433323436)
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 2]),0.37170414232138)
        #check if sign of eigenvectors is consistent
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 0]))
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 1]))
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 2]))
        #test irreps
        self.assertEqual(self.Contributions._IRLabels[-1],'B2')

        #TODO: test NAC
        # self.assertAlmostEqual(self.Contributions2._frequencies[47], 3490.6434922723,places=1)
        # #test calculation of eigenvectors
        # self.assertAlmostEqual(abs(self.Contributions2._EigFormat[15, 47, 0]),0.00084433323436)
        # self.assertAlmostEqual(abs(self.Contributions2._EigFormat[15, 47, 1]),0.00084433323436)
        # self.assertAlmostEqual(abs(self.Contributions2._EigFormat[15, 47, 2]),0.37170414232138)
        # #check if sign of eigenvectors is consistent
        # self.assertEqual(np.sign(self.Contributions2._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 0]))
        # self.assertEqual(np.sign(self.Contributions2._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 1]))
        # self.assertEqual(np.sign(self.Contributions2._EigFormat[14, 47, 2]),np.sign(self.Contributions._EigFormat[15, 47, 2]))
        # #test irreps
        # self.assertEqual(self.Contributions2._IRLabels[-1],'B2')



        #TODO: set masses externally [e.g., use D mass]

        #TODO: start from FORCE constants instead


if __name__ == '__main__':
    unittest.main()


