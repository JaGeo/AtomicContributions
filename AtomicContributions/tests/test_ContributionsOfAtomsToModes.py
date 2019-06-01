from AtomicContributions.ContributionsOfAtomsToModes import AtomicContributionsCalculator

import unittest
import numpy as np
import os

path_here = os.path.dirname(__file__)


class AtomicContributionToModesTest(unittest.TestCase):
    def setUp(self):
        self.Contributions = AtomicContributionsCalculator(PoscarName=os.path.join(path_here, 'POSCAR'),
                                                           ForceConstants=False,
                                                           ForceFileName=os.path.join(path_here, 'FORCE_SETS'),
                                                           supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],
                                                           primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.Contributions_masses = AtomicContributionsCalculator(PoscarName=os.path.join(path_here, 'POSCAR'),
                                                                  ForceConstants=False,
                                                                  ForceFileName=os.path.join(path_here, 'FORCE_SETS'),
                                                                  supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],
                                                                  primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                                  masses=[12.010700, 12.010700, 15.999400, 15.999400,
                                                                          14.006700, 14.006700, 14.006700, 14.006700, 2,
                                                                          2,
                                                                          2, 2, 2, 2, 2, 2])

        self.Contributions2 = AtomicContributionsCalculator(PoscarName=os.path.join(path_here, 'POSCAR.NaCl'),
                                                            ForceConstants=False,
                                                            ForceFileName=os.path.join(path_here, 'FORCE_SETS.NaCl'),
                                                            supercell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]], nac=True,
                                                            BornFileName=os.path.join(path_here, 'BORN.NaCl'),
                                                            primitive=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
        self.ContributionsFC = AtomicContributionsCalculator(PoscarName=os.path.join(path_here, 'POSCAR_Methanol'),
                                                             ForceConstants=True,
                                                             ForceFileName=os.path.join(path_here,
                                                                                        'FORCE_CONSTANTS_Methanol'),
                                                             supercell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], nac=False)

    def test_attributes(self):
        # test calculation of frequencies
        self.assertAlmostEqual(self.Contributions._frequencies[47], 3490.6434922723, places=1)
        # test calculation of eigenvectors
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 0]), 0.00084433323436)
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 1]), 0.00084433323436)
        self.assertAlmostEqual(abs(self.Contributions._EigFormat[15, 47, 2]), 0.37170414232138)
        # check if sign of eigenvectors is consistent!!
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions._EigFormat[15, 47, 0]))
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions._EigFormat[15, 47, 1]))
        self.assertEqual(np.sign(self.Contributions._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions._EigFormat[15, 47, 2]))
        # test irreps
        self.assertEqual(self.Contributions._IRLabels[-1], 'B2')
        # test contributions
        sum_contribution = 0.0
        for atom in range(0, 16):
            sum_contribution += self.Contributions._PercentageAtom[47, atom]
        self.assertAlmostEqual(sum_contribution, 1.0)

        # TODO: test NAC
        self.assertAlmostEqual(self.Contributions2._frequencies[-1], 153.7212069157, places=2)

        # TODO: set masses externally [e.g., use D mass]
        self.assertAlmostEqual(self.Contributions_masses._frequencies[47], 2598.2875793589, places=1)
        # test calculation of eigenvectors
        self.assertAlmostEqual(abs(self.Contributions_masses._EigFormat[15, 47, 0]), 0.00378948635566)
        self.assertAlmostEqual(abs(self.Contributions_masses._EigFormat[15, 47, 1]), 0.00378948635566)
        self.assertAlmostEqual(abs(self.Contributions_masses._EigFormat[15, 47, 2]), 0.33223420830758)
        # check if sign of eigenvectors is consistent
        self.assertEqual(np.sign(self.Contributions_masses._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions_masses._EigFormat[15, 47, 0]))
        self.assertEqual(np.sign(self.Contributions_masses._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions_masses._EigFormat[15, 47, 1]))
        self.assertEqual(np.sign(self.Contributions_masses._EigFormat[14, 47, 2]),
                         np.sign(self.Contributions_masses._EigFormat[15, 47, 2]))
        # test irreps
        self.assertEqual(self.Contributions._IRLabels[-1], 'B2')

        # start from FORCE constants instead
        self.assertAlmostEqual(self.ContributionsFC._frequencies[-1], 3741.4132865293, places=1)


if __name__ == '__main__':
    unittest.main()
