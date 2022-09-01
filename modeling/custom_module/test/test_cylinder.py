import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm
import sympy as sm


def point_lin_dist(p1, p2, p3):
    # line goes through p1 and p2
    p1 = sm.Point(*p1)
    p2 = sm.Point(*p2)
    l = sm.Line(p1, p2)
    return float(l.distance(sm.Point(*p3)))


m = IMP.Model()
particle_coordinates = [
    (0, 0, 0),
    (0, 0, 1),
    (0, 3, 0),
    (5, 0, 0),
    (10, 8, 7),
    (3, 6, 9),
    (-4, -2.5, -4)
]
centers = [
    (0, 0, 0),
    (0, 0, 1),
    (0, 3, 0),
    (5, 0, 0),
    (3, 9, 6),
    (-4, -2.5, -4)
]
radii = [1, 2, 2.5, 3]
kappas = [1, 4]


def setup_p(name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


particles = [setup_p(i, particle_coordinates[i]) for i in range(len(particle_coordinates))]
# check the non-general formulation (parallel to axis)
for axis in [0, 1, 2]:
    for k in kappas:
        for c in centers:
            for r in radii:
                c1, c2 = [c[i] for i in range(len(c)) if i != axis]
                res = desm.CylinderLocalizationRestraint(particles, k, axis, c1, c2, r)
                val = res.unprotected_evaluate(None)
                new_c = list(c)
                new_c[axis] += 1
                answer_val = [point_lin_dist(c, new_c, x) for x in particle_coordinates]
                answer_val = [(x - r) for x in answer_val if x > r]
                answer_val = sum([k * (x ** 2) for x in answer_val])
                match = abs(val - answer_val) < 1e-7
                print(f'Special -> {match}')


# check the general formulation (parallel to axis)
for c1 in centers:
    for k in kappas:
        for c2 in centers:
            if c1 == c2:
                continue
            # Only works if the axis is not already parallel to z axis
            if (c1[0] == c2[0]) and (c1[1] == c2[1]):
                continue
            for r in radii:
                res = desm.CylinderLocalizationRestraint(particles, k, c1[0], c1[1], c1[2], c2[0], c2[1], c2[2], r)
                val = res.unprotected_evaluate(None)
                answer_val = [point_lin_dist(c1, c2, x) for x in particle_coordinates]
                answer_val = [(x - r) for x in answer_val if x > r]
                answer_val = sum([k * (x ** 2) for x in answer_val])
                match = abs(val - answer_val) < 1e-7
                print(f'General -> {match}')