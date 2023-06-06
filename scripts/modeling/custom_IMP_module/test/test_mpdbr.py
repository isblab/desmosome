import IMP
import IMP.core
import IMP.algebra
import IMP.pmi
import IMP.pmi.restraints
import IMP.container
from scipy.spatial.distance import cdist
import numpy as np
import IMP.test


# TODO: Add x0 != 0 test case
def setup_p(m, name, vec, r):
    p = m.add_particle(str(name))
    p = IMP.core.XYZR.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    p.set_radius(r)
    return p


class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist1, plist2, x0=0, kappa=0, label=None, weight=1):
        name = 'MinimumPairDistanceBindingRestraint%1%'
        super(MinimumPairDistanceBindingRestraint, self).__init__(model, name=name, label=label, weight=weight)
        l1 = IMP.container.ListSingletonContainer(model)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(model)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added MPDBR on particles", len(plist1), ":", len(plist2), "at x0:kappa", str(x0), ":",
              str(kappa))


class Tests(IMP.test.TestCase):
    """ TestMinimum Pair Distance Binding Restraint Restraint"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (0, 0, 0),
            (0, 0, 1),
            (0, 3, 0),
            (5, 0, 0),
            (10, 8, 7),
            (3, 6, 9),
            (-4, -2.5, -4)
        ]
        self.kp = 2
        self.radii = [1, 0.5, 3, 2, 0.75, 1.1, 0.1]

        dist_matrix = cdist(np.array(self.particle_coordinates[:3]),
                            np.array(self.particle_coordinates[3:]),
                            metric='euclidean')
        radii_sums = np.array(self.radii[:3])[:, np.newaxis] + np.array(self.radii[3:])[np.newaxis, :]
        dist_matrix -= radii_sums
        self.answer1 = (dist_matrix.flatten()[np.argmin(np.abs(dist_matrix.flatten()))] ** 2) * self.kp / 2

        dist_matrix = cdist(np.array(self.particle_coordinates[::2]),
                            np.array(self.particle_coordinates[1::2]),
                            metric='euclidean')
        radii_sums = np.array(self.radii[::2])[:, np.newaxis] + np.array(self.radii[1::2])[np.newaxis, :]
        dist_matrix -= radii_sums
        self.answer2 = (dist_matrix.flatten()[np.argmin(np.abs(dist_matrix.flatten()))] ** 2) * self.kp / 2

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i], self.radii[i])
                          for i in range(len(self.particle_coordinates))]

    def test_mpdbr(self):
        """Test Minimum Pair Distance Binding Restraint"""
        res1 = MinimumPairDistanceBindingRestraint(self.m, self.particles[:3], self.particles[3:], 0, self.kp)
        res2 = MinimumPairDistanceBindingRestraint(self.m, self.particles[::2], self.particles[1::2], 1, self.kp)
        self.assertAlmostEqual(self.answer1, res1.evaluate(), 7)
        self.assertAlmostEqual(self.answer2, res2.evaluate(), 7)


if __name__ == '__main__':
    IMP.test.main()
