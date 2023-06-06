import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm
import math
import IMP.test


def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


class Tests(IMP.test.TestCase):
    """ Test Single-axis Min-Gaussian Restraint"""

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
        self.means = [0, 1, -1.5, 2, 4]
        self.sigmas = [0.5, 5, 10, 3]
        self.particles = [setup_p(self.m, i, self.particle_coordinates[i])
                          for i in range(len(self.particle_coordinates))]

    def test_gaussian(self):
        """Test Single-axis Min-Gaussian restraint with multiple sigma/kappa/means"""
        for axis in [0, 1, 2]:
            for mn in self.means:
                for s in self.sigmas:
                    res = desm.SingleAxisMinGaussianRestraint(self.particles, axis, mn, s)
                    val = res.unprotected_evaluate(None)
                    answer_val = min([abs(p[axis] - mn) for p in self.particle_coordinates])
                    answer_val = math.log(2 * math.pi * (s ** 2)) / 2 + (answer_val ** 2) / (2 * (s ** 2))
                    self.assertAlmostEqual(val, answer_val, 7)


if __name__ == '__main__':
    IMP.test.main()
