import RMF
import sys
import warnings
import re
from collections import defaultdict
import time
import os
import IMP
import IMP.pmi
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.topology
import IMP.pmi.tools
import IMP.pmi.dof
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import random
import math
import string
from fixed_shuffle import shuffle_configuration

random.seed(12345)

try:
    import IMP.desmosome
except ModuleNotFoundError:
    warnings.warn("The Script needs the Desmosome Module!", UserWarning)


# FUNCTIONS AND WRAPPERS -----------------------------------------------


class SingleAxisMinGaussianRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, mean, sigma, axis, label, weight=1):
        particles = plist
        name = 'SingleAxisMinGaussianRestraint%1%'
        super(SingleAxisMinGaussianRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.SingleAxisMinGaussianRestraint(particles, axis, mean, sigma)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added SAMGR on particles", len(plist), "at mean:sigma:axis", str(mean), ":", str(sigma), ":",
              str(axis))


class CylinderLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, axis, center1, center2, radius, kappa, label, weight=1):
        particles = plist
        name = 'CylinderLocalizationRestraint%1%'
        super(CylinderLocalizationRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.CylinderLocalizationRestraint(particles, kappa, axis, center1, center2, radius)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added CLR on particles", len(plist), "at center1:center2:r:kappa", str(center1), ":",
              str(center2), ":", str(radius), ":", str(kappa))


class AxialLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, axis, axisUpper, axisLower, kappa, label, weight=1):
        particles = plist
        name = 'AxialLocalizationRestraint%1%'
        super(AxialLocalizationRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.AxialLocalizationRestraint(particles, axis, axisUpper, axisLower, kappa)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added ALR on particles", len(plist), "at axisU:axisL:kappa", str(axisUpper), ":",
              str(axisLower), ":", str(kappa))


def chain_id_gen():  # To sequentially generate 52 alphabets to use as Chain IDs
    for i in (list(string.ascii_uppercase)):
        yield i
    for i in (list(string.ascii_lowercase)):
        yield i


def param_parse():
    m = [re.match("-[A-Z_]+=.+", i) for i in sys.argv[1:]]
    ind = 1
    given_params = []
    for mm in m:
        if mm:
            val = sys.argv[ind]
            given_params.append(val[1:])
        ind += 1
    return given_params


def set_polarity_flip_y(molecule, cval_is_greater=True):
    leaves = IMP.atom.get_leaves(molecule)
    leaves = [x for x in leaves if
              IMP.core.RigidMember.get_is_setup(x.get_particle()) or IMP.core.NonRigidMember.get_is_setup(
                  x.get_particle())]  # get only the structured parts
    nval = IMP.core.XYZ(leaves[0].get_particle()).get_coordinate(1)
    cval = IMP.core.XYZ(leaves[-1].get_particle()).get_coordinate(1)
    if cval_is_greater:  # The C-region is supposed to be away from the plasma membrane
        check = (nval > cval)
    else:  # The N-region is supposed to be away from the plasma membrane
        check = (cval > nval)
    if check:  # If a flip is necessary to restore the "correct" polarity
        rb, temp = IMP.pmi.tools.get_rbs_and_beads(molecule)
        assert len(rb) == 1, "Multiple Rigid Bodies passed for polarity flip"
        rb = rb[0]
        identity_rotation = IMP.algebra.get_identity_rotation_3d()
        identity_translation = IMP.algebra.Vector3D(0, 0, 0)
        current_position = IMP.core.XYZ(rb).get_coordinates()
        random_axis = IMP.algebra.Vector3D(random.random(), 0, random.random())
        magnitude = math.sqrt(random_axis[0] ** 2 + random_axis[2] ** 2)
        random_axis /= magnitude
        flip_rotation = IMP.algebra.Rotation3D(0, random_axis[0], random_axis[1], random_axis[2])
        IMP.core.transform(rb, IMP.algebra.Transformation3D(identity_rotation, -current_position))
        IMP.core.transform(rb, IMP.algebra.Transformation3D(flip_rotation, identity_translation))
        IMP.core.transform(rb, IMP.algebra.Transformation3D(identity_rotation, current_position))
        return True  # A flip was performed
    return False  # No flip needed


# FUNCTIONS AND WRAPPERS -----------------------------------------------

# FLAGS AND PARAMETERS -----------------------------------------------
post_vars = None  # To read this variable in the computation below
pre_vars = set(locals().keys())  # All the locally defined variable names

RB_COLLISION = False  # Whether to attempt to avoid collisions while randomizing
DP_NUMBER = 2
PG_NUMBER = 2
PG_FLAG_NC = True  # To add PG Beads to EM density or not
DP_FLAG_NC = True  # To add DP Beads to EM density or not
PG_MT = 2  # PG Rigid Body max trans (Includes the structured DC region)
DP_MT = 1  # DP Rigid Body max trans
BEAD_MT = 4  # Beads max trans
DP_S_MT = 0.4  # DP Super Rigid Body max trans
PG_S_MT = 0.7  # PG Super Rigid Body max trans (Includes the structured DC region)
CONN_WEIGHT = 10  # Weight of Connectivity Restraint
CONN_SCALE_PG = 2.68  # Connectivity Restraint Scale factor for PG
CONN_SCALE_DP = 1.89  # Connectivity Restraint Scale factor for DP
CONN_RES = 20  # Resolution to apply the connectivity restraint at
RB_RES = [30]  # Resolutions for the representation of the structured parts (Rigid Bodies) for all molecules
PG_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in PG
DP_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in DP
EVR_STRUCTURE_WEIGHT = 10  # Weight of EVR Restraint for structured regions
EVR_BEAD_WEIGHT = 10  # Weight of EVR Restraint for Beads (optionally separate from above in implementation)
EVR_KAPPA = 1  # EVR Restraint Kappa factor
EM_WEIGHT_PGDP = 100  # Weight of EM Restraint for PGDP layer (optionally separate from above in implementation)
EM_SLOPE_PGDP = 1e-8  # Slope parameter for the Gaussian EM Restraint for the PG-DP molecules
SAMGR_RES = 20  # Resolution of the representation to apply the SAMGR Restraint on
SAMGR_NPG = 10  # Weight of SAMGR Restraint for N-PG
SAMGR_CPG = 10  # Weight of SAMGR Restraint for C-PG
SAMGR_NDP = 10  # Weight of SAMGR Restraint for N-DP
CYL_WEIGHT = 10  # Weight of Cylinder Localization Restraint
AXIAL_WEIGHT = 10  # Weight of Axial Localization Restraint
REX_MAX_TEMP = 2.5  # Maximum temperature for Replica Exchange
OUTPUT_PATH = "output"  # For all the RMFs and STAT files
AXIAL_FLIP_PGDP = True  # Whether to flip the PG/DP to ensure the correct relative positioning of C/N terminals
CYLINDER_RADIUS = 135  # Radius of the cylinder for the Cylindrical Localization restraint
NFRAMES = 7500
PG_COLOR = "#E6194B"
DP_COLOR = "#FFE119"

post_vars = set(locals().keys())  # All the locally defined variable names
parameter_variables = post_vars - pre_vars  # Only the parameters defined above
params_to_override = param_parse()
for p in params_to_override:  # Only override the parameters defined above (guards against random typos too)
    if p.split("=")[0] not in parameter_variables:
        print("WARNING: Option " + p + " overrides a non-run parameter. Value is neglected for now")
    else:
        print("PARSER: Executing option", p)
        exec(p.replace('#', '\"'))

dens = "pgdp_separate_gmm_20.txt"
if not os.path.isdir(OUTPUT_PATH):  # make the output folder if none exist
    os.mkdir(OUTPUT_PATH)
# FLAGS AND PARAMETERS -----------------------------------------------

# REPRESENTATION -----------------------------------------------
sequences = dict()  # load the FASTA sequences
for i in [x for x in os.listdir("data") if ".fasta" in x]:
    sequences[i.split("_uniprot")[0]] = IMP.pmi.topology.Sequences("data/" + i)

start_time = time.time()

mdl = IMP.Model()
syst = IMP.pmi.topology.System(mdl)
st = syst.create_state()

pg_molecules = []
dp_molecules = []

chainGen = chain_id_gen()
# Create the molecules and their copies (Currently all chains are named sequentially
pg_molecules.append(st.create_molecule("PG", sequences["pg"]["pg"], chain_id=chainGen.__next__()))
for i in range(PG_NUMBER - 1):
    pg_molecules.append(pg_molecules[0].create_copy(chain_id=chainGen.__next__()))
dp_molecules.append(st.create_molecule("DP", sequences["dp"]["dp"], chain_id=chainGen.__next__()))
for i in range(DP_NUMBER - 1):
    dp_molecules.append(dp_molecules[0].create_copy(chain_id=chainGen.__next__()))

# Add structures based on the PDBs for each of the copies
structures = defaultdict(list)  # Initializes to an empty list for a new key
for i, p in enumerate(dp_molecules):
    x = p.add_structure("data/3r6n.pdb", chain_id="A", res_range=(178, 584))
    structures['dp'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/dp_" + str(i),
                         density_voxel_size=6, color=DP_COLOR)
    p.add_representation(p[:584] - x, DP_FB_RES, setup_particles_as_densities=DP_FLAG_NC, color=DP_COLOR)
for i, p in enumerate(pg_molecules):
    x = p.add_structure("data/pg-dsc1-pg-dsc1_modified.pdb", chain_id="A", res_range=(126, 673))
    structures['pg'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/pg_" + str(i),
                         density_voxel_size=6, color=PG_COLOR)
    p.add_representation(p[:] - x, PG_FB_RES, setup_particles_as_densities=PG_FLAG_NC, color=PG_COLOR)

root_hierarchy = syst.build()
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
# REPRESENTATION -----------------------------------------------

# DEGREES OF FREEDOM -----------------------------------------------
rigid_bodies = []
super_rigid_bodies = []
flexible_beads = []
for i in range(DP_NUMBER):
    s = dp_molecules[i][177:]
    s2 = dp_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s, max_trans=DP_MT, max_rot=0.1, nonrigid_max_trans=BEAD_MT, nonrigid_parts=s & s2,
                              name="dp_structure_" + str(i))
    rigid_bodies.append(x)
    x = dof.create_flexible_beads(dp_molecules[i][:] - s, max_trans=BEAD_MT)
    flexible_beads.append(x)
    x = dof.create_super_rigid_body(dp_molecules[i][:], max_trans=DP_S_MT, max_rot=0.1, name="dp_srb_" + str(i))
    super_rigid_bodies.append(x)
for i in range(PG_NUMBER):
    s = pg_molecules[i][125:673]
    s2 = pg_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s, max_trans=PG_MT, max_rot=0.1, nonrigid_max_trans=BEAD_MT, nonrigid_parts=s & s2,
                              name="pg_dsc1_structure_" + str(i))
    rigid_bodies.append(x)
    x = dof.create_flexible_beads(pg_molecules[i][:] - s, max_trans=BEAD_MT)
    flexible_beads.append(x)
    x = dof.create_super_rigid_body(pg_molecules[i][:], max_trans=PG_S_MT, max_rot=0.1,
                                    name="pg_dsc1_srb_" + str(i))
    super_rigid_bodies.append(x)
# DEGREES OF FREEDOM -----------------------------------------------

print('REPRESENTATION AND DOF LOAD TIME: ', round(time.time() - start_time, 2))
if '--topology_only' in sys.argv:
    quit(0)

# RESTRAINTS -------------------------------------------------------------
restraint_list = []

# connectivity restraint
scale_mapping = {'PG': CONN_SCALE_PG, 'DP': CONN_SCALE_DP}
for m in root_hierarchy.get_children()[0].get_children():  # get the molecules in state 0
    scale = scale_mapping[m.get_name()]
    cr_kwargs = {'scale': scale, 'disorderedlength': False, 'upperharmonic': True, 'resolution': CONN_RES,
                 'label': m.get_name() + '_'}
    connectivity_restraint = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m, **cr_kwargs)
    connectivity_restraint.set_weight(CONN_WEIGHT)
    connectivity_restraint.add_to_model()
    restraint_list.append(connectivity_restraint)

# excluded volume restraint
evr_kwargs = {'kappa': EVR_KAPPA, 'resolution': 1000}
excluded_volume_restraint = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=[root_hierarchy],
                                                                                    **evr_kwargs)
# TODO: Separate this into beads and structures
excluded_volume_restraint.set_weight(EVR_STRUCTURE_WEIGHT)
excluded_volume_restraint.add_to_model()
restraint_list.append(excluded_volume_restraint)

# SHUFFLING -----------------------------------------------
# Output a RMF frame to see the validity of the representation and the connectivity/excluded volume restraints
molecules = root_hierarchy.get_children()[0].get_children()
span = 1.1 * CYLINDER_RADIUS * math.sqrt(2)
bb_pg_dp = ((141 - span / 2, 245, 141 - span / 2), (141 + span / 2, 345, 141 + span / 2))

temp = [x for x in molecules if ('PG' == x.get_name()) or ('DP' == x.get_name())]
try:
    shuffle_configuration(temp, bounding_box=bb_pg_dp, verbose=True, niterations=500, avoidcollision_rb=RB_COLLISION)
except ValueError:
    print("SHUFFLING: Failed once. Restarting with a lower distance cutoff")
    shuffle_configuration(temp, bounding_box=bb_pg_dp, verbose=True, niterations=500, cutoff=3,
                          avoidcollision_rb=RB_COLLISION)
print("SHUFFLING: Shuffling Done!")

f = RMF.create_rmf_file(OUTPUT_PATH + '/shuffling_and_optimization_' + str(os.getpid()) + '.rmf')
IMP.rmf.add_hierarchy(f, root_hierarchy)
IMP.rmf.save_frame(f)
if AXIAL_FLIP_PGDP:
    pgs = [x for x in molecules if x.get_name() == 'PG']
    dps = [x for x in molecules if x.get_name() == 'DP']
    for i in range(PG_NUMBER):  # Set the polarity of the PGs and DPs to be appropriate
        flipped = set_polarity_flip_y(pgs[i], False)
        if flipped:
            print("SHUFFLING: Flipped PG molecule ", i)
    for i in range(DP_NUMBER):
        flipped = set_polarity_flip_y(dps[i], True)
        if flipped:
            print("SHUFFLING: Flipped DP molecule ", i)
    IMP.rmf.save_frame(f)
print("SHUFFLING: Optimizing flexible Beads.")
dof.optimize_flexible_beads(nsteps=1000, temperature=2)
IMP.rmf.save_frame(f)
del f
# SHUFFLING -----------------------------------------------

print('SHUFFLING TIME: ', round(time.time() - start_time, 2))
if '--shuffle_only' in sys.argv:
    quit(0)

# EM Restraint
density_particles = IMP.atom.Selection(root_hierarchy, molecules=['PG', 'DP'],
                                       representation_type=IMP.atom.DENSITIES).get_selected_particles()
em_kwargs = {'target_fn': './data/gmm/' + dens, 'slope': EM_SLOPE_PGDP, 'scale_target_to_mass': True,
             'weight': EM_WEIGHT_PGDP}
em_restraint = IMP.pmi.restraints.em.GaussianEMRestraint(density_particles, **em_kwargs)
em_restraint.set_label("fullPGDP")
em_restraint.add_to_model()
restraint_list.append(em_restraint)

# ImmunoGold restraint
# PG
for i in range(PG_NUMBER):
    selection_tuple = (1, 106, 'PG', i, None)
    plist = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, SAMGR_RES)
    selection_tuple2 = (666, 738, 'PG', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple2, SAMGR_RES)
    yGaussian = SingleAxisMinGaussianRestraintC(mdl, plist, 115 + 229, 4.5, 1, 'YGaussianNTermPG' + str(i), SAMGR_NPG)
    yGaussian.add_to_model()
    restraint_list.append(yGaussian)
    yGaussian = SingleAxisMinGaussianRestraintC(mdl, plist2, 115 + 108, 9, 1, 'YGaussianCTermPG' + str(i), SAMGR_CPG)
    yGaussian.add_to_model()
    restraint_list.append(yGaussian)
# DP
for i in range(DP_NUMBER):
    selection_tuple = (1, 189, 'DP', i, None)
    plist = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, SAMGR_RES)
    yGaussian = SingleAxisMinGaussianRestraintC(mdl, plist, 115 + 103, 9.8, 1, 'YGaussianNTermDP' + str(i), SAMGR_NDP)
    yGaussian.add_to_model()
    restraint_list.append(yGaussian)

# Localization restraints
all_particles = IMP.atom.Selection(root_hierarchy, resolution=1000,
                                   representation_type=IMP.atom.BALLS).get_selected_particles()
cylRes = CylinderLocalizationRestraintC(mdl, all_particles, 1, 132, 132, CYLINDER_RADIUS, 5,
                                        "CylinderLocalizationMain", CYL_WEIGHT)
cylRes.add_to_model()
restraint_list.append(cylRes)

axialRes = AxialLocalizationRestraintC(mdl, all_particles, 1, 430, 180, 5, "AxialLocalizationMain", AXIAL_WEIGHT)
axialRes.add_to_model()
restraint_list.append(axialRes)
# RESTRAINTS -----------------------------------------------

print('SETUP TIME: ', round(time.time() - start_time, 2))
if '--setup_only' in sys.argv:
    quit(0)

rex = IMP.pmi.macros.ReplicaExchange0(mdl,
                                      root_hier=root_hierarchy,
                                      monte_carlo_sample_objects=dof.get_movers(),
                                      global_output_directory=OUTPUT_PATH,
                                      output_objects=restraint_list,
                                      monte_carlo_steps=10,
                                      replica_exchange_maximum_temperature=REX_MAX_TEMP,
                                      replica_exchange_minimum_temperature=1,
                                      number_of_best_scoring_models=0,
                                      number_of_frames=NFRAMES)
rex.execute_macro()

print('TOTAL TIME: ', round(time.time() - start_time, 2))
quit(0)
