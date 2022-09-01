import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm
import math

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
means = [0, 1, -1.5, 2, 4]
sigmas = [0.5, 5, 10, 3]


def setup_p(name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


particles = [setup_p(i, particle_coordinates[i]) for i in range(len(particle_coordinates))]
for axis in [0, 1, 2]:
    for mn in means:
        for s in sigmas:
            res = desm.SingleAxisMinGaussianRestraint(particles, axis, mn, s)
            val = res.unprotected_evaluate(None)
            answer_val = min([abs(p[axis] - mn) for p in particle_coordinates])
            answer_val = math.log(2 * math.pi * (s ** 2)) / 2 + (answer_val ** 2) / (2 * (s ** 2))
            match = abs(val - answer_val) < 1e-7
            print(f'axis: {axis} + mean: {mn} + sigma: {s} -> {match}')