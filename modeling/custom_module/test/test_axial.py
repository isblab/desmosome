import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm
import IMP.test


def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


class Tests(IMP.test.TestCase):
    """ Test Axial Localization Restraint"""

    def setUp(self):
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
        self.upper_list = [1, 2, 2, 4, 6, -2, -3, 0]
        self.lower_list = [-1, 0, -4, 3.5, 1, -6, -3.5, -2]
        self.kappas = [1, 4]
        self.answers = {0: [110, 90, 74, 130.25, 44, 230, 296.25, 138],
                        1: [80.25, 59.25, 53, 93, 19.25, 201, 265.25, 109.25],
                        2: [109, 90, 74, 133.25, 38, 223, 287.25, 135]}

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i])
                          for i in range(len(self.particle_coordinates))]

    def test_axial(self):
        """Test Axial Restrain with multiple axes/kappa"""
        for axis in [0, 1, 2]:
            for k in self.kappas:
                count = 0
                for i, j in zip(self.upper_list, self.lower_list):
                    res = desm.AxialLocalizationRestraint(self.particles, axis, i, j, k)
                    val = res.unprotected_evaluate(None)
                    self.assertAlmostEqual(val / k - self.answers[axis][count], 0, 1e-7)
                    count += 1


if __name__ == '__main__':
    IMP.test.main()
