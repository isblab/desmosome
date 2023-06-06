import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm
import IMP.test
import sympy as sm


def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


def point_lin_dist(p1, p2, p3):
    # line goes through p1 and p2
    p1 = sm.Point(*p1)
    p2 = sm.Point(*p2)
    l = sm.Line(p1, p2)
    return float(l.distance(sm.Point(*p3)))


class Tests(IMP.test.TestCase):
    """ Test Cylinder Localization Restraint"""

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
        self.centers = [
            (0, 0, 0),
            (0, 0, 1),
            (0, 3, 0),
            (5, 0, 0),
            (3, 9, 6),
            (-4, -2.5, -4)
        ]
        self.radii = [1, 2, 2.5, 3]
        self.kappas = [1, 4]

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def test_cylinder_special(self):
        """Test Cylinder Localization Restraint special version (axis-parallel)"""

        # check the non-general formulation (parallel to axis)
        for axis in [0, 1, 2]:
            for k in self.kappas:
                for c in self.centers:
                    for r in self.radii:
                        c1, c2 = [c[i] for i in range(len(c)) if i != axis]
                        res = desm.CylinderLocalizationRestraint(self.particles, k, axis, c1, c2, r)
                        val = res.unprotected_evaluate(None)
                        new_c = list(c)
                        new_c[axis] += 1
                        answer_val = [point_lin_dist(c, new_c, x) for x in self.particle_coordinates]
                        answer_val = [(x - r) for x in answer_val if x > r]
                        answer_val = sum([k * (x ** 2) for x in answer_val])
                        self.assertAlmostEqual(val, answer_val, 7)

    def test_cylinder_general(self):
        """Test Cylinder Localization Restraint general version"""
        # check the general formulation (parallel to axis)
        for c1 in self.centers:
            for k in self.kappas:
                for c2 in self.centers:
                    if c1 == c2:
                        continue
                    # Only works if the axis is not already parallel to z axis
                    if (c1[0] == c2[0]) and (c1[1] == c2[1]):
                        continue
                    for r in self.radii:
                        res = desm.CylinderLocalizationRestraint(self.particles, k, c1[0], c1[1],
                                                                 c1[2], c2[0], c2[1], c2[2], r)
                        val = res.unprotected_evaluate(None)
                        answer_val = [point_lin_dist(c1, c2, x) for x in self.particle_coordinates]
                        answer_val = [(x - r) for x in answer_val if x > r]
                        answer_val = sum([k * (x ** 2) for x in answer_val])
                        self.assertAlmostEqual(val, answer_val, 7)


if __name__ == '__main__':
    IMP.test.main()

