import IMP
import IMP.core
import IMP.pmi
import IMP.pmi.tools
import IMP.core
import math


def shuffle_configuration(objects,
                          max_translation=300., max_rotation=2.0 * math.pi,
                          avoidcollision_rb=True, avoidcollision_fb=False,
                          cutoff=10.0, niterations=100,
                          bounding_box=None,
                          excluded_rigid_bodies=[],
                          hierarchies_excluded_from_collision=[],
                          hierarchies_included_in_collision=[],
                          verbose=False,
                          return_debug=False):
    """Shuffle particles. Used to restart the optimization.
    The configuration of the system is initialized by placing each
    rigid body and each bead randomly in a box. If `bounding_box` is
    specified, the particles are placed inside this box; otherwise, each
    particle is displaced by up to max_translation angstroms, and randomly
    rotated. Effort is made to place particles far enough from each other to
    prevent any steric clashes.
    @param objects Can be one of the following inputs:
               IMP Hierarchy, PMI System/State/Molecule/TempResidue, or
               a list/set of them
    @param max_translation Max translation (rbs and flexible beads)
    @param max_rotation Max rotation (rbs only)
    @param avoidcollision_rb check if the particle/rigid body was
           placed close to another particle; uses the optional
           arguments cutoff and niterations
    @param avoidcollision_fb Advanced. Generally you want this False because
           it's hard to shuffle beads.
    @param cutoff Distance less than this is a collision
    @param niterations How many times to try avoiding collision
    @param bounding_box Only shuffle particles within this box.
           Defined by ((x1,y1,z1),(x2,y2,z2)).
    @param excluded_rigid_bodies Don't shuffle these rigid body objects
    @param hierarchies_excluded_from_collision Don't count collision
           with these bodies
    @param hierarchies_included_in_collision Hierarchies that are not
           shuffled, but should be included in collision calculation
           (for fixed regions)
    @param verbose Give more output
    @note Best to only call this function after you've set up degrees
          of freedom
    For debugging purposes, returns: <shuffled indexes>,
    <collision avoided indexes>
    """

    # checking input
    hierarchies = IMP.pmi.tools.input_adaptor(objects,
                                              pmi_resolution='all',
                                              flatten=True)
    rigid_bodies, flexible_beads = IMP.pmi.tools.get_rbs_and_beads(hierarchies)
    if len(rigid_bodies) > 0:
        mdl = rigid_bodies[0].get_model()
    elif len(flexible_beads) > 0:
        mdl = flexible_beads[0].get_model()
    else:
        raise Exception("Could not find any particles in the hierarchy")
    if len(rigid_bodies) == 0:
        print("shuffle_configuration: rigid bodies were not initialized")

    # gather all particles
    gcpf = IMP.core.GridClosePairsFinder()
    gcpf.set_distance(cutoff)

    # Add particles from excluded hierarchies to excluded list
    collision_excluded_hierarchies = IMP.pmi.tools.input_adaptor(
        hierarchies_excluded_from_collision, pmi_resolution='all',
        flatten=True)

    collision_included_hierarchies = IMP.pmi.tools.input_adaptor(
        hierarchies_included_in_collision, pmi_resolution='all', flatten=True)

    collision_excluded_idxs = set(
        leaf.get_particle().get_index()
        for h in collision_excluded_hierarchies
        for leaf in IMP.core.get_leaves(h))

    collision_included_idxs = set(
        leaf.get_particle().get_index()
        for h in collision_included_hierarchies
        for leaf in IMP.core.get_leaves(h))

    # Excluded collision with Gaussians
    all_idxs = []  # expand to representations?
    for p in IMP.pmi.tools.get_all_leaves(hierarchies):
        if IMP.core.XYZ.get_is_setup(p):
            all_idxs.append(p.get_particle_index())
        if IMP.core.Gaussian.get_is_setup(p):
            collision_excluded_idxs.add(p.get_particle_index())

    if bounding_box is not None:
        ((x1, y1, z1), (x2, y2, z2)) = bounding_box
        ub = IMP.algebra.Vector3D(x1, y1, z1)
        lb = IMP.algebra.Vector3D(x2, y2, z2)
        bb = IMP.algebra.BoundingBox3D(ub, lb)

    all_idxs = set(all_idxs) | collision_included_idxs
    all_idxs = all_idxs - collision_excluded_idxs
    debug = []
    print('shuffling', len(rigid_bodies), 'rigid bodies')
    for rb in rigid_bodies:
        if rb not in excluded_rigid_bodies:
            # gather particles to avoid with this transform
            if avoidcollision_rb:
                rb_idxs = set(rb.get_member_particle_indexes()) - \
                          collision_excluded_idxs
                other_idxs = all_idxs - rb_idxs

            # iterate, trying to avoid collisions
            niter = 0
            while niter < niterations:
                rbxyz = (rb.get_x(), rb.get_y(), rb.get_z())

                # local transform
                if bounding_box:
                    translation = IMP.algebra.get_random_vector_in(bb)
                    # First move to origin
                    transformation_orig = IMP.algebra.Transformation3D(
                        IMP.algebra.get_identity_rotation_3d(),
                        -IMP.core.XYZ(rb).get_coordinates())
                    IMP.core.transform(rb, transformation_orig)
                    rotation = IMP.algebra.get_random_rotation_3d()
                    transformation = IMP.algebra.Transformation3D(rotation,
                                                                  translation)

                else:
                    transformation = \
                        IMP.algebra.get_random_local_transformation(
                            rbxyz, max_translation, max_rotation)

                debug.append([rb, other_idxs if avoidcollision_rb else set()])
                IMP.core.transform(rb, transformation)

                # check collisions
                if avoidcollision_rb and other_idxs:
                    mdl.update()
                    npairs = len(gcpf.get_close_pairs(mdl,
                                                      list(other_idxs),
                                                      list(rb_idxs)))
                    if npairs == 0:
                        break
                    else:
                        niter += 1
                        if verbose:
                            print("shuffle_configuration: rigid body placed "
                                  "close to other %d particles, trying "
                                  "again..." % npairs)
                            print("shuffle_configuration: rigid body name: "
                                  + rb.get_name())
                        if niter == niterations:
                            raise ValueError(
                                "tried the maximum number of iterations to "
                                "avoid collisions, increase the distance "
                                "cutoff")
                else:
                    break

    print('shuffling', len(flexible_beads), 'flexible beads')
    for fb in flexible_beads:
        # gather particles to avoid
        if avoidcollision_fb:
            fb_idxs = set(IMP.get_indexes([fb]))
            other_idxs = all_idxs - fb_idxs
            if not other_idxs:
                continue

        # iterate, trying to avoid collisions
        niter = 0
        while niter < niterations:
            if bounding_box:
                translation = IMP.algebra.get_random_vector_in(bb)
                transformation = IMP.algebra.Transformation3D(translation)
            else:
                fbxyz = IMP.core.XYZ(fb).get_coordinates()
                transformation = IMP.algebra.get_random_local_transformation(
                    fbxyz, max_translation, max_rotation)

            # For gaussians, treat this fb as an rb
            if IMP.core.NonRigidMember.get_is_setup(fb):
                memb = IMP.core.NonRigidMember(fb)
                xyz = memb.get_internal_coordinates()
                if bounding_box:
                    # 'translation' is the new desired position in global
                    # coordinates; we need to convert that to internal
                    # coordinates first using the rigid body's ref frame
                    rf = memb.get_rigid_body().get_reference_frame()
                    glob_to_int = rf.get_transformation_from()
                    memb.set_internal_coordinates(
                            glob_to_int.get_transformed(translation))
                else:
                    xyz_transformed = transformation.get_transformed(xyz)
                    memb.set_internal_coordinates(xyz_transformed)
                debug.append([xyz, other_idxs if avoidcollision_fb else set()])
            else:
                d = IMP.core.XYZ(fb)
                if bounding_box:
                    # Translate to origin first
                    if IMP.core.RigidBody.get_is_setup(fb.get_particle()):
                        IMP.core.transform(
                            IMP.core.RigidBody(fb.get_particle()),
                            -d.get_coordinates())
                    else:
                        IMP.core.transform(d, -d.get_coordinates())
                    d = IMP.core.XYZ(fb)
                debug.append([d, other_idxs if avoidcollision_fb else set()])
                if IMP.core.RigidBody.get_is_setup(fb.get_particle()):
                    IMP.core.transform(
                        IMP.core.RigidBody(fb.get_particle()), transformation)
                else:
                    IMP.core.transform(d, transformation)

            if avoidcollision_fb:
                mdl.update()
                npairs = len(gcpf.get_close_pairs(mdl,
                                                  list(other_idxs),
                                                  list(fb_idxs)))

                if npairs == 0:
                    break
                else:
                    niter += 1
                    print("shuffle_configuration: floppy body placed close "
                          "to other %d particles, trying again..." % npairs)
                    if niter == niterations:
                        raise ValueError(
                            "tried the maximum number of iterations to avoid "
                            "collisions, increase the distance cutoff")
            else:
                break
    if return_debug:
        return debug