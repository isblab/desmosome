import IMP
import IMP.core
import IMP.algebra
import IMP.pmi
import IMP.pmi.restraints
import IMP.container
from scipy.spatial.distance import cdist
import numpy as np

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
kp = 2
radii = [1, 0.5, 3, 2, 0.75, 1.1, 0.1]

dist_matrix = cdist(np.array(particle_coordinates[:3]), np.array(particle_coordinates[3:]), metric='euclidean')
radii_sums = np.array(radii[:3])[:, np.newaxis] + np.array(radii[3:])[np.newaxis, :]
dist_matrix -= radii_sums
answer1 = (dist_matrix.flatten()[np.argmin(np.abs(dist_matrix.flatten()))] ** 2) * kp / 2

dist_matrix = cdist(np.array(particle_coordinates[::2]), np.array(particle_coordinates[1::2]), metric='euclidean')
radii_sums = np.array(radii[::2])[:, np.newaxis] + np.array(radii[1::2])[np.newaxis, :]
dist_matrix -= radii_sums
answer2 = (dist_matrix.flatten()[np.argmin(np.abs(dist_matrix.flatten()))] ** 2) * kp / 2


def setup_p(name, vec, r):
    p = m.add_particle(str(name))
    p = IMP.core.XYZR.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    p.set_radius(r)
    return p


particles = [setup_p(i, particle_coordinates[i], radii[i]) for i in range(len(particle_coordinates))]


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


res1 = MinimumPairDistanceBindingRestraint(m, particles[:3], particles[3:], 0, kp)
res2 = MinimumPairDistanceBindingRestraint(m, particles[::2], particles[1::2], 0, kp)
print(f'Test 1: {np.abs(answer1 - res1.evaluate()) < 1e-7}')
print(f'Test 2: {np.abs(answer2 - res2.evaluate()) < 1e-7}')