import IMP
import IMP.core
import IMP.algebra
import IMP.desmosome as desm

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
upper_list = [1, 2, 2, 4, 6, -2, -3, 0]
lower_list = [-1, 0, -4, 3.5, 1, -6, -3.5, -2]
kappas = [1, 4]
answers = {0: [110, 90, 74, 130.25, 44, 230, 296.25, 138],
           1: [80.25, 59.25, 53, 93, 19.25, 201, 265.25, 109.25],
           2: [109, 90, 74, 133.25, 38, 223, 287.25, 135]}


def setup_p(name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p


particles = [setup_p(i, particle_coordinates[i]) for i in range(len(particle_coordinates))]
for axis in [0, 1, 2]:
    for k in kappas:
        count = 0
        for i, j in zip(upper_list, lower_list):
            res = desm.AxialLocalizationRestraint(particles, axis, i, j, k)
            val = res.unprotected_evaluate(None)
            match = abs(val / k - answers[axis][count]) < 1e-7
            print(f'axis: {axis} + kappa: {k} + count: {count} -> {match}')
            count += 1