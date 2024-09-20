import RMF
import sys
import warnings
from collections import defaultdict
import time
import os
import IMP
import IMP.pmi
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.samplers
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
from IMP.pmi.tools import shuffle_configuration
import pickle

# Imports needed to use ProtocolOutput
import IMP.pmi.mmcif
import ihm
import ihm.location
import ihm.model


random.seed(12345)

try:
    import IMP.desmosome
except ModuleNotFoundError:
    warnings.warn("The Script needs the Desmosome Module!", UserWarning)


# TODO: Double check the flipping code and the methods to get only the structure (use returned rigid bodies ^ directly?)

# FUNCTIONS AND WRAPPERS -----------------------------------------------
# wrapper for the Binding Restraint
class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist1, plist2, x0=0, kappa=1, label=None, weight=1):
        name = 'MinimumPairDistanceBindingRestraint%1%'
        super(MinimumPairDistanceBindingRestraint, self).__init__(model, name=name, label=label, weight=weight)
        l1 = IMP.container.ListSingletonContainer(mdl)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(mdl)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added MPDBR on particles", len(plist1), ":", len(plist2), "at x0:kappa", str(x0), ":",
              str(kappa))


# wrapper for the SAMGR defined in custom c++ code
class SingleAxisMinGaussianRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, mean, sigma, axis, label, weight=1):
        particles = plist
        name = 'SingleAxisMinGaussianRestraint%1%'
        super(SingleAxisMinGaussianRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.SingleAxisMinGaussianRestraint(particles, axis, mean, sigma)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added SAMGR on particles", len(plist), "at mean:sigma:axis", str(mean), ":", str(sigma), ":",
              str(axis))


# wrapper for the Cylinder Localization Restraint defined in custom c++ code
class CylinderLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, axis, center1, center2, radius, kappa, label, weight=1):
        particles = plist
        name = 'CylinderLocalizationRestraint%1%'
        super(CylinderLocalizationRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.CylinderLocalizationRestraint(particles, kappa, axis, center1, center2, radius)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added CLR on particles", len(plist), "at center1:center2:r:kappa", str(center1), ":",
              str(center2), ":", str(radius), ":", str(kappa))


# To create a BallMover which only has 2 axes' FloatKeys setup to prevent moving in the axis_to_drop direction
def add_axially_restricted_ball(dof_obj, selection_obj, max_trans):
    assert isinstance(dof_obj, IMP.pmi.dof.DegreesOfFreedom)
    # Currently hard-coding to drop the movement along y-axis
    float_key_array = [IMP.FloatKey('x'), IMP.FloatKey('z')]
    # The below code goes along the layout of the code for pmi.add_flexible_beads with appropriate modifications
    hiers = IMP.pmi.tools.input_adaptor(selection_obj, "all", flatten=True)
    assert len(hiers) == 1, 'Give only a single bead'
    p = hiers[0].get_particle()
    IMP.core.XYZ(p).set_coordinates_are_optimized(True)
    if IMP.core.RigidMember.get_is_setup(hiers[0]) or IMP.core.NonRigidMember.get_is_setup(hiers[0]):
        raise Exception("Cannot create flexible beads from members of rigid body")
    dof_obj.flexible_beads.append(hiers[0])
    fbmv = IMP.core.BallMover(p.get_model(), p, float_key_array, max_trans)
    fb_movers = [fbmv]
    fbmv.set_was_used(True)
    dof_obj.fb_movers.append(fbmv)
    dof_obj.movers_particles_map[fbmv] = IMP.atom.get_leaves(hiers[0])
    dof_obj.movers_xyz_map[fbmv] = IMP.atom.get_leaves(hiers[0])
    dof_obj.movers += fb_movers
    return fb_movers


def chain_id_gen():  # To sequentially generate 52 alphabets to use as Chain IDs
    for i in (list(string.ascii_uppercase)):
        yield i
    for i in (list(string.ascii_lowercase)):
        yield i


# to flip DP and PG (structured regions only) if their N-C polarity is reversed
def set_polarity_flip_y(molecule, cval_is_greater=True):
    # this assumes that the structured region does not loop back. For generality,
    # the max and min y-coordinate should be taken instead of the nval and cval below
    leaves = IMP.atom.get_leaves(molecule)  # the lowest resolution in the hierarchy
    leaves = [x for x in leaves if
              IMP.core.RigidMember.get_is_setup(x.get_particle()) or IMP.core.NonRigidMember.get_is_setup(
                  x.get_particle())]  # get only the structured parts
    # assuming the leaves array is populated in an ordered fashion
    # get the coordinates of the first and the last bead (the N-C termini of the structured part)
    nval = IMP.core.XYZ(leaves[0].get_particle()).get_coordinate(1)
    cval = IMP.core.XYZ(leaves[-1].get_particle()).get_coordinate(1)
    if cval_is_greater:  # the C-region is supposed to be away from the plasma membrane
        check = (nval > cval)
    else:  # the N-region is supposed to be away from the plasma membrane
        check = (cval > nval)
    if check:  # tf a flip is necessary to restore the "correct" polarity
        rb, temp = IMP.pmi.tools.get_rbs_and_beads(molecule)
        assert len(rb) == 1, "Multiple Rigid Bodies passed for polarity flip"
        rb = rb[0]
        identity_rotation = IMP.algebra.get_identity_rotation_3d()
        identity_translation = IMP.algebra.Vector3D(0, 0, 0)
        current_position = IMP.core.XYZ(rb).get_coordinates()
        # pick any axis in the x-z plane (such that rotation along this of 180 flips the y-polarity)
        random_axis = IMP.algebra.Vector3D(random.random(), 0, random.random())
        magnitude = math.sqrt(random_axis[0] ** 2 + random_axis[2] ** 2)
        random_axis /= magnitude  # unit vector
        # create a representative quaternion for rotation
        # rotation of angle t around the axis (ux, uy, uz) gives a quaternion
        # (cos(t/2), (ux, uy, uz) * sin(t/2)). Here, rotating by 180 degrees
        flip_rotation = IMP.algebra.Rotation3D(0, random_axis[0], random_axis[1], random_axis[2])
        # shift to origin -> rotate -> shift back to original position
        IMP.core.transform(rb, IMP.algebra.Transformation3D(identity_rotation, -current_position))
        IMP.core.transform(rb, IMP.algebra.Transformation3D(flip_rotation, identity_translation))
        IMP.core.transform(rb, IMP.algebra.Transformation3D(identity_rotation, current_position))
        return True  # flip was performed
    return False  # no flip needed


# to calculate the connectivity-restraint scale for disordered proteins
def flory_self_avoiding(num_residues, num_segments, bead_radius, base_conn_dist=3.6):
    rg = 1.92 * (num_residues ** 0.6)  # based on Kohn PNAS
    a = math.sqrt((rg ** 2) * 176 / 25 / (num_segments ** 1.2))  # from Flory
    # alternatively, assume <(r_i - r_j)^2> = b^2 * |i - j|^(2v) with v = 6/5
    total_length = a * num_segments
    scale = ((total_length - 2 * bead_radius) / num_segments - 2 * bead_radius) / base_conn_dist
    return int(scale)  # truncating to the nearest smaller integer
# FUNCTIONS AND WRAPPERS -----------------------------------------------


# FLAGS AND PARAMETERS -----------------------------------------------
# three early-stopping flags --topology_only, --shuffle_only, --setup_only

PG_FLAG_NC = True  # To add PG Beads to EM density or not
DP_FLAG_NC = True  # To add DP Beads to EM density or not
PKP_FLAG_NC = False  # To add PKP beads to EM density or not
PKP_MT = 3.25  # PKP Rigid Body max trans
PG_MT = 0.4  # PG Rigid Body max trans (Includes the structured DC region)
DP_MT = 0.4  # DP Rigid Body max trans
BEAD_DC_MT = 7  # Beads max trans for DSC/DSG
BEAD_PKP_MT = 8  # Beads max trans for PKP
BEAD_PG_MT = 7  # Beads max trans for PG
BEAD_DP_MT = 5  # Beads max trans for DP
PG_S_MT = 0.3  # PG Super Rigid Body max trans (Includes the structured DC region)
RB_RES = [30]  # Resolutions for the representation of the structured parts (Rigid Bodies) for all molecules
DC_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in DC
PKP_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in PKP
PG_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in PG
DP_FB_RES = 20  # Resolution for the representation of the unstructured flexible bead in DP
PKP_MR = 0.12  # PKP Rigid Body max rot
PG_MR = 0.065  # PG Rigid Body max rot (Includes the structured DC region)
DP_MR = 0.1  # DP Rigid Body max rot
PG_S_MR = 0.025  # PG Super Rigid Body max rot (Includes the structured DC region)
CONN_WEIGHT = 100  # Weight of Connectivity Restraint
CONN_SCALE = 'auto'  # The scale parameter for the connectivity restraint
CONN_RES = 20  # Resolution to apply the connectivity restraint at
EVR_STRUCTURE_WEIGHT = 10  # Weight of EVR Restraint for structured regions
EVR_BEAD_WEIGHT = 10  # Weight of EVR Restraint for Beads (optionally separate from above in implementation)
EVR_KAPPA = 1  # EVR Restraint Kappa factor
EM_WEIGHT_PKP = 100  # Weight of EM Restraint for PKP layer
EM_WEIGHT_PGDP = 100  # Weight of EM Restraint for PGDP layer (optionally separate from above in implementation)
EM_SLOPE_PKP = 1e-8  # Slope parameter for the Gaussian EM Restraint for the PKP molecules
EM_SLOPE_PGDP = 1e-8  # Slope parameter for the Gaussian EM Restraint for the PG-DP molecules
SAMGR_RES = 20  # Resolution of the representation to apply the SAMGR Restraint on
SAMGR_NPKP = 10  # Weight of SAMGR Restraint for N-PKP
SAMGR_CPKP = 10  # Weight of SAMGR Restraint for C-PKP
SAMGR_NPG = 10  # Weight of SAMGR Restraint for N-PG
SAMGR_CPG = 10  # Weight of SAMGR Restraint for C-PG
SAMGR_NDP = 10  # Weight of SAMGR Restraint for N-DP
MPDBR_WEIGHT = 5  # Weight of MPDBR Restraint
COLOC_FACTOR = 1  # Factor to multiply the MPDBR weight to weigh co-localization experimental evidence
OVERLAY_FACTOR = 3  # Factor to multiply the MPDBR weight to weigh overlay assay experimental evidence
Y2H_FACTOR = 2  # Factor to multiply the MPDBR weight to weigh Y2H experimental evidence
COIP_FACTOR = 3  # Factor to multiply the MPDBR weight to weigh co-IP experimental evidence
MPDBR_RES = 20  # Resolution of the representation to apply the MPDBR Restraint on
CYL_WEIGHT = 100  # Weight of Cylinder Localization Restraint
PKP_BLOB_SIZE = 5  # Half the length of the diagonal of the blob-wise bounding box for shuffling
REX_MAX_TEMP = 2.5  # Maximum temperature for Replica Exchange
OUTPUT_PATH = sys.argv[1]  # For all the RMFs and STAT files
AXIAL_FLIP_PGDP = True  # Whether to flip the PG/DP to ensure the correct relative positioning of C/N terminals
CYLINDER_RADIUS = 150  # Radius of the cylinder for the Cylindrical Localization restraint
RB_COLLISION = False  # Whether to attempt to avoid collisions while randomizing
PKP_NUMBER = 4  # How many actual PKP molecules to keep
GPKP_NUMBER = 3  # How many "ghost" PKP molecules to keep
DP_NUMBER = 4  # How many DP molecules to keep
PG_NUMBER = 4  # How many PG molecules to keep -> same number of (DSCs + DSGs)

NFRAMES = 2 #Note only 2 frames for mmcif
C_PKP_COLOR = '#95d293'
STRUCTURE_PKP_COLOR = '#4daf4a'  # Includes the color for GPKPs
N_PKP_COLOR = '#377e35'
DP_COLOR = '#e31a1c'
DSC_COLOR = '#1f78b4'
DSG_COLOR = '#6a3d9a'
N_PG_COLOR = '#994d00'
STRUCT_PG_COLOR = '#ff7f00'
C_PG_COLOR = '#ffb366'

# use layer wise densities
dens = ["pkp_separate_gmm_10.txt", "pgdp_separate_gmm_20.txt"]

if not os.path.isdir(OUTPUT_PATH):  # make the output folder if none exist
    os.mkdir(OUTPUT_PATH)
assert PKP_NUMBER == DP_NUMBER, "Can't setup the PKP-DP binding restraint for unequal PKP/DP"
assert PG_NUMBER == DP_NUMBER, "Can't setup the PG-DP binding restraint for unequal PG/DP"
assert GPKP_NUMBER > 0, "The script is hard-coded to assume at least 1 GPKP"
# FLAGS AND PARAMETERS -----------------------------------------------

# REPRESENTATION -----------------------------------------------
sequences = dict()  # load the FASTA sequences
for i in [x for x in os.listdir("data") if ".fasta" in x]:
    sequences[i.split("_uniprot")[0]] = IMP.pmi.topology.Sequences("data/" + i)

start_time = time.time()

mdl = IMP.Model()
syst = IMP.pmi.topology.System(mdl)
st = syst.create_state()

# Add deposition information
po = IMP.pmi.mmcif.ProtocolOutput()
syst.add_protocol_output(po)
po.system.title = "Integrative structure of the epithelial desmosomal outer plaque"
# po.system.citations.append(ihm.Citation.from_pubmed_id(000000)) #TODO

pkp_molecules = []
gpkp_molecules = []
pg_molecules = []
dp_molecules = []
dsc_molecules = []
dsg_molecules = []

chainGen = chain_id_gen()
# Create the molecules and their copies (Currently all chains are named sequentially)
# TODO: hard code the numbers
pkp_molecules.append(st.create_molecule("PKP1a", sequences["pkp1a"]["pkp1a"], chain_id=chainGen.__next__()))
for i in range(PKP_NUMBER - 1):
    pkp_molecules.append(pkp_molecules[0].create_copy(chain_id=chainGen.__next__()))
gpkp_molecules.append(st.create_molecule("GPKP1a", sequences["pkp1a"]["pkp1a"], chain_id=chainGen.__next__()))
for i in range(GPKP_NUMBER - 1):
    gpkp_molecules.append(gpkp_molecules[0].create_copy(chain_id=chainGen.__next__()))
pg_molecules.append(st.create_molecule("PG", sequences["pg"]["pg"], chain_id=chainGen.__next__()))
for i in range(PG_NUMBER - 1):
    pg_molecules.append(pg_molecules[0].create_copy(chain_id=chainGen.__next__()))
dp_molecules.append(st.create_molecule("DP", sequences["dp"]["dp"], chain_id=chainGen.__next__()))
for i in range(DP_NUMBER - 1):
    dp_molecules.append(dp_molecules[0].create_copy(chain_id=chainGen.__next__()))
dsc_molecules.append(st.create_molecule("DSC1", sequences["dsc1"]["dsc1"], chain_id=chainGen.__next__()))
for i in range(PG_NUMBER - (PG_NUMBER // 2) - 1):
    dsc_molecules.append(dsc_molecules[0].create_copy(chain_id=chainGen.__next__()))
dsg_molecules.append(st.create_molecule("DSG1", sequences["dsg1"]["dsg1"], chain_id=chainGen.__next__()))
for i in range(PG_NUMBER // 2 - 1):
    dsg_molecules.append(dsg_molecules[0].create_copy(chain_id=chainGen.__next__()))

# Add structures based on the PDBs for each of the copies
structures = defaultdict(list)
for i, p in enumerate(pkp_molecules):  # The "real" PKPs
    x = p.add_structure("data/1xm9.pdb", chain_id="A", res_range=(244, 700))
    structures['pkp1a'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/pkp1a_" + str(i),
                         density_voxel_size=6, color=STRUCTURE_PKP_COLOR)
    p.add_representation(p[243:700] - x, PKP_FB_RES, setup_particles_as_densities=True, color=STRUCTURE_PKP_COLOR)
    p.add_representation(p[:243], PKP_FB_RES, setup_particles_as_densities=PKP_FLAG_NC, color=N_PKP_COLOR)
    p.add_representation(p[700:], PKP_FB_RES, setup_particles_as_densities=PKP_FLAG_NC, color=C_PKP_COLOR)
for i, p in enumerate(gpkp_molecules):  # The "ghost" PKPs
    x = p.add_structure("data/docked_" + str(i + 1) + "_gpkp.pdb", chain_id="A", res_range=(244, 700))
    structures['gpkp1a'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/gpkp1a_" + str(i),
                         density_voxel_size=6, color=STRUCTURE_PKP_COLOR)
for i, p in enumerate(dp_molecules):
    x = p.add_structure("data/3r6n.pdb", chain_id="A", res_range=(178, 584))
    structures['dp'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/dp_" + str(i),
                         density_voxel_size=6, color=DP_COLOR)
    p.add_representation(p[:584] - x, DP_FB_RES, setup_particles_as_densities=DP_FLAG_NC, color=DP_COLOR)
dsc_number = len(dsc_molecules)
for i, p in enumerate(pg_molecules[:dsc_number]):
    x = p.add_structure("data/pg-dsc1-pg-dsc1_modified.pdb", chain_id="A", res_range=(126, 673))
    structures['pg'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/pg_" + str(i),
                         density_voxel_size=6, color=STRUCT_PG_COLOR)
    p.add_representation(p[670:] - x, PG_FB_RES, setup_particles_as_densities=PG_FLAG_NC, color=C_PG_COLOR)
    p.add_representation(p[:128] - x, PG_FB_RES, setup_particles_as_densities=PG_FLAG_NC, color=N_PG_COLOR)
    p.add_representation(p[128:670] - x, PG_FB_RES, setup_particles_as_densities=True, color=STRUCT_PG_COLOR)
for i, p in enumerate(pg_molecules[dsc_number:]):
    x = p.add_structure("data/pg-dsg1-pg-dsg1_modified.pdb", chain_id="A", res_range=(126, 673))
    structures['pg'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10,
                         density_prefix="data/gmm/pg_" + str(i + dsc_number),
                         density_voxel_size=6, color=STRUCT_PG_COLOR)
    p.add_representation(p[670:] - x, PG_FB_RES, setup_particles_as_densities=PG_FLAG_NC, color=C_PG_COLOR)
    p.add_representation(p[:128] - x, PG_FB_RES, setup_particles_as_densities=PG_FLAG_NC, color=N_PG_COLOR)
    p.add_representation(p[128:670] - x, PG_FB_RES, setup_particles_as_densities=True, color=STRUCT_PG_COLOR)
for i, p in enumerate(dsc_molecules):
    x = p.add_structure("data/pg-dsc1-pg-dsc1_modified.pdb", chain_id="C", res_range=(834, 894))
    structures['dsc'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/dsc1_" + str(i),
                         density_voxel_size=6, color=DSC_COLOR)
    p.add_representation(p[714:] - x, DC_FB_RES, setup_particles_as_densities=False, color=DSC_COLOR)
for i, p in enumerate(dsg_molecules):
    x = p.add_structure("data/pg-dsg1-pg-dsg1_modified.pdb", chain_id="C", res_range=(698, 765))
    structures['dsg'].append(x)
    p.add_representation(x, RB_RES, density_residues_per_component=10, density_prefix="data/gmm/dsg1_" + str(i),
                         density_voxel_size=6, color=DSG_COLOR)
    p.add_representation(p[569:842] - x, DC_FB_RES, setup_particles_as_densities=False, color=DSG_COLOR)

root_hierarchy = syst.build()


dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
# REPRESENTATION -----------------------------------------------

# DEGREES OF FREEDOM -----------------------------------------------
rigid_bodies = []
super_rigid_bodies = []
flexible_beads = []
movers_in_order = []  # used while analyzing the runs for identifying the protein-specific bead movers
for i in range(PKP_NUMBER):
    s = pkp_molecules[i][243:700]
    s2 = pkp_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s, max_trans=PKP_MT, max_rot=PKP_MR, nonrigid_max_trans=BEAD_PKP_MT, nonrigid_parts=s & s2,
                              name="pkp1a_structure_" + str(i))
    rigid_bodies.append(x)
    movers_in_order += list(zip(x[0], ['pkp1a'] * len(x[0])))  # TODO: Comment
    x = dof.create_flexible_beads(pkp_molecules[i][:] - s, max_trans=BEAD_PKP_MT)
    movers_in_order += list(zip(x, ['pkp1a'] * len(x)))
    flexible_beads.append(x)
for i in range(GPKP_NUMBER):
    s = gpkp_molecules[i][243:700]
    s2 = gpkp_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s - s2, max_trans=0, max_rot=0, name="gpkp1a_structure_" + str(i))
    movers_in_order += list(zip(x[0], ['gpkp1a'] * len(x[0])))
    rigid_bodies.append(x)
for i in range(DP_NUMBER):
    s = dp_molecules[i][177:]
    s2 = dp_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s, max_trans=DP_MT, max_rot=DP_MR, nonrigid_max_trans=BEAD_DP_MT, nonrigid_parts=s & s2,
                              name="dp_structure_" + str(i))
    movers_in_order += list(zip(x[0], ['dp'] * len(x[0])))
    rigid_bodies.append(x)
    x = dof.create_flexible_beads(dp_molecules[i][:] - s, max_trans=BEAD_DP_MT)
    movers_in_order += list(zip(x, ['dp'] * len(x)))
    flexible_beads.append(x)
for i in range(len(dsc_molecules)):
    s = pg_molecules[i][125:673]
    s2 = pg_molecules[i].get_non_atomic_residues()
    s3 = dsc_molecules[i][833:]
    s3_nar = dsc_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s | s3, max_trans=PG_MT, max_rot=PG_MR, nonrigid_max_trans=BEAD_PG_MT,
                              nonrigid_parts=(s & s2) | (s3 & s3_nar), name="pg_dsc1_structure_" + str(i))
    movers_in_order += list(zip(x[0], ['pg_dsc1'] * len(x[0])))
    rigid_bodies.append(x)
    x = dof.create_flexible_beads(pg_molecules[i][:] - s, max_trans=BEAD_PG_MT)
    movers_in_order += list(zip(x, ['pg'] * len(x)))
    flexible_beads.append(x)
    x = dof.create_flexible_beads(dsc_molecules[i][714 + DC_FB_RES:] - s3, max_trans=BEAD_DC_MT)
    movers_in_order += list(zip(x, ['dsc1'] * len(x)))
    flexible_beads.append(x)
    x = add_axially_restricted_ball(dof, dsc_molecules[i][714:714 + DC_FB_RES], BEAD_DC_MT)
    movers_in_order += list(zip(x, ['dsc1'] * len(x)))
    flexible_beads.append(x)
    x = dof.create_super_rigid_body(s | dsc_molecules[i][714 + DC_FB_RES:], max_trans=PG_S_MT, max_rot=PG_S_MR,
                                    name="pg_dsc1_srb_" + str(i))
    movers_in_order += list(zip(x, ['pg_dsc1'] * len(x)))
    super_rigid_bodies.append(x)
for i in range(len(dsg_molecules)):
    s = pg_molecules[i + len(dsc_molecules)][125:673]
    s2 = pg_molecules[i + len(dsc_molecules)].get_non_atomic_residues()
    s3 = dsg_molecules[i][697:765]
    s3_nar = dsg_molecules[i].get_non_atomic_residues()
    x = dof.create_rigid_body(s | s3, max_trans=PG_MT, max_rot=PG_MR, nonrigid_max_trans=BEAD_PG_MT,
                              nonrigid_parts=(s & s2) | (s3 & s3_nar), name="pg_dsg1_structure_" + str(i))
    movers_in_order += list(zip(x[0], ['pg_dsg1'] * len(x[0])))
    rigid_bodies.append(x)
    x = dof.create_flexible_beads(pg_molecules[i + len(dsc_molecules)][:] - s, max_trans=BEAD_PG_MT)
    flexible_beads.append(x)
    movers_in_order += list(zip(x, ['pg'] * len(x)))
    x = dof.create_flexible_beads(dsg_molecules[i][569 + DC_FB_RES:842] - s3, max_trans=BEAD_DC_MT)
    flexible_beads.append(x)
    movers_in_order += list(zip(x, ['dsg1'] * len(x)))
    x = add_axially_restricted_ball(dof, dsg_molecules[i][569:569 + DC_FB_RES], BEAD_DC_MT)
    flexible_beads.append(x)
    movers_in_order += list(zip(x, ['dsg1'] * len(x)))
    x = dof.create_super_rigid_body(s | dsg_molecules[i][569 + DC_FB_RES:842], max_trans=PG_S_MT, max_rot=PG_S_MR,
                                    name="pg_dsg1_srb_" + str(i))
    movers_in_order += list(zip(x, ['pg_dsg1'] * len(x)))
    super_rigid_bodies.append(x)

movers_in_order = [(str(x[0]).replace("\"", ""), str(x[1])) for x in movers_in_order]
# DEGREES OF FREEDOM -----------------------------------------------

print('REPRESENTATION AND DOF LOAD TIME: ', round(time.time() - start_time, 2))
if '--topology_only' in sys.argv:
    f = open('mover_list_list_object', 'wb')
    pickle.dump(movers_in_order, f)
    f.close()
    quit(0)

# RESTRAINTS -------------------------------------------------------------
# any random bead can be used of appropriate size
first_particle = IMP.atom.Selection(root_hierarchy, molecule='DSC1', residue_indexes=[715]).get_selected_particles()
radius = IMP.core.XYZR(first_particle[0]).get_radius()
print("BEAD_RADIUS:", radius)
if CONN_SCALE == 'auto':
    CONN_SCALE_DC = flory_self_avoiding(125, 125 // DC_FB_RES, radius)  # average of DSG/DSC num_residues
    CONN_SCALE_PKP = flory_self_avoiding(244, 244 // PKP_FB_RES, radius)
    # weighted average based on psipred predicted disordered regions
    CONN_SCALE_PG = (106 * flory_self_avoiding(126, 126 // DP_FB_RES, radius) + 20) / 126
    CONN_SCALE_DP = (53 * flory_self_avoiding(178, 178 // PG_FB_RES, radius) + 125) / 178
else:
    assert len(CONN_SCALE) == 4, 'Improper CONN_SCALE (only 4-element lists are acceptable)'
    CONN_SCALE_DC, CONN_SCALE_PKP, CONN_SCALE_PG, CONN_SCALE_DP = CONN_SCALE
restraint_list = []

# connectivity restraint
scale_mapping = {'PKP1a': CONN_SCALE_PKP, 'PG': CONN_SCALE_PG, 'DP': CONN_SCALE_DP, 'DSC1': CONN_SCALE_DC,
                 'DSG1': CONN_SCALE_DC, 'GPKP1a': CONN_SCALE_PKP}
count = 0  # just to differentiate the different molecules' connectivity restraints
for m in root_hierarchy.get_children()[0].get_children():  # get the molecules in state 0
    scale = scale_mapping[m.get_name()]
    cr_kwargs = {'scale': scale, 'disorderedlength': False, 'upperharmonic': True, 'resolution': CONN_RES,
                 'label': m.get_name() + '_' + str(count) + '_'}
    connectivity_restraint = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m, **cr_kwargs)
    connectivity_restraint.set_weight(CONN_WEIGHT)
    connectivity_restraint.add_to_model()
    restraint_list.append(connectivity_restraint)
    count += 1

# excluded volume restraint
evr_kwargs = {'kappa': EVR_KAPPA, 'resolution': 1000}
excluded_volume_restraint = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=[root_hierarchy],
                                                                                    **evr_kwargs)
excluded_volume_restraint.set_weight(EVR_STRUCTURE_WEIGHT)
excluded_volume_restraint.add_to_model()
restraint_list.append(excluded_volume_restraint)

# SHUFFLING -----------------------------------------------
# Output a RMF frame to see the validity of the representation and the connectivity/excluded volume restraints
# the pkp and gpkp blobs represent the blob-centers of the respective EM blobs
# 115A along the y-axis is the location of the upper surface of the pm
molecules = root_hierarchy.get_children()[0].get_children()
bb_pkp_layer = ((15, 130, 15), (255, 220, 255))
bb_pkp_blobs = [(72.967, 162.808, 142.779), (188.050, 156.253, 191.862), (137.378, 163.380, 143.345),
                (176.426, 158.955, 95.430)]
# bb_gpkp_blobs is unused (written here only to include all the blobs)
bb_gpkp_blobs = [(107.886, 168.625, 191.741), (223.253, 154.809, 134.812), (110.415, 164.555, 83.391)]
# side of the maximal square inside the cylinder
span = CYLINDER_RADIUS * math.sqrt(2)
# get the bounding boxes
bb_pg_dp = ((141 - span / 2, 245, 141 - span / 2), (141 + span / 2, 345, 141 + span / 2))
bb_dc = ((141 - span / 2, 130, 141 - span / 2), (141 + span / 2, 360, 141 + span / 2))
bb_anchor = ((141 - span / 2, 115 + radius, 141 - span / 2), (141 + span / 2, 115 + radius + 2, 141 + span / 2))
# bb_anchor ensures the lower edge of the bead is within 2A of the pm
temp = [x for x in molecules if ('DSC1' == x.get_name()) or ('DSG1' == x.get_name())]
shuffle_configuration(temp, bounding_box=bb_dc, avoidcollision_rb=False)
temp = [x for x in molecules if ('PG' == x.get_name()) or ('DP' == x.get_name())]
try:  # often fails since the layer is crowded if RB_COLLISION is True
    # if the RB_COLLISION is False, it does not reach the except clause usually
    shuffle_configuration(temp, bounding_box=bb_pg_dp, verbose=True, niterations=5, avoidcollision_rb=RB_COLLISION) #Note:  decreased the number of shuffle iterations
except ValueError:
    print("SHUFFLING: Failed once. Restarting with a lower distance cutoff")
    shuffle_configuration(temp, bounding_box=bb_pg_dp, verbose=True, niterations=5, cutoff=3,
                          avoidcollision_rb=RB_COLLISION) #Note:  decreased the number of shuffle iterations

temp = [x for x in molecules if 'PKP1a' == x.get_name()]
temp_g = [x for x in molecules if 'GPKP1a' == x.get_name()]  # no need to shuffle GPKP
shuffle_configuration(temp, bounding_box=bb_pkp_layer, avoidcollision_rb=False)
# shuffling to randomize the orientation since it is easier than to randomly rotate it explicitly below
for i, j in zip(temp, bb_pkp_blobs):  # for each rigid body, place its center in one of the PKP blobs
    rbs, fbs = IMP.pmi.tools.get_rbs_and_beads(i)
    assert len(rbs) == 1, 'More than one RB passed'
    # The center of the molecule is randomized in a small box around the EM-blob center
    # The following code is taken (with small modifications) from pmi shuffle
    k1 = (j[0] - PKP_BLOB_SIZE, j[1] - PKP_BLOB_SIZE, j[2] - PKP_BLOB_SIZE)
    k2 = (j[0] + PKP_BLOB_SIZE, j[1] + PKP_BLOB_SIZE, j[2] + PKP_BLOB_SIZE)
    ((x1, y1, z1), (x2, y2, z2)) = k1, k2  # unpacking the coordinates
    ub = IMP.algebra.Vector3D(x1, y1, z1)
    lb = IMP.algebra.Vector3D(x2, y2, z2)
    bb = IMP.algebra.BoundingBox3D(ub, lb)
    translation = IMP.algebra.get_random_vector_in(bb)  # random location of the center in the bb
    translation -= IMP.core.XYZ(rbs[0]).get_coordinates()
    transform = IMP.algebra.Transformation3D(IMP.algebra.get_identity_rotation_3d(), translation)
    IMP.core.transform(rbs[0], transform)

# only the first Bead in the DSC and DSG are anchored
anchored_residues = [IMP.atom.get_leaves(x)[0] for x in molecules if 'DSC1' == x.get_name()]
anchored_residues += [IMP.atom.get_leaves(x)[0] for x in molecules if 'DSG1' == x.get_name()]
shuffle_configuration(anchored_residues, bounding_box=bb_anchor, avoidcollision_rb=False)
print("SHUFFLING: Shuffling Done!")

# output some frames to a debug-RMF to check for shuffling, axial-flip and bead-optimization
f = RMF.create_rmf_file(OUTPUT_PATH + '/shuffling_and_optimization_' + str(os.getpid()) + '.rmf')
IMP.rmf.add_hierarchy(f, root_hierarchy)
IMP.rmf.save_frame(f)
if AXIAL_FLIP_PGDP:  # TODO: Remove this
    pgs = [x for x in molecules if x.get_name() == 'PG']
    dps = [x for x in molecules if x.get_name() == 'DP']
    for i in range(PG_NUMBER):  # Set the polarity of the PGs and DPs to be appropriate
        flipped = set_polarity_flip_y(pgs[i], False)
        if flipped:
            print("SHUFFLING: Flipped PG molecule ", i)
        flipped = set_polarity_flip_y(dps[i], True)
        if flipped:
            print("SHUFFLING: Flipped DP molecule ", i)
    IMP.rmf.save_frame(f)
print("SHUFFLING: Optimizing flexible Beads.")
# dof.optimize_flexible_beads(nsteps=1000, temperature=1)
# setting up a custom optimizer instead of dof.optimize_flexible_beads to avoid the anchored residues from getting
# a separate unrestricted mover
mc = IMP.pmi.samplers.MonteCarlo(dof.model, dof.get_floppy_body_movers(), 1)
assert len(dof.get_floppy_body_movers()) == len(dof.get_flexible_beads()), 'Number of movers != number of fbs'
print('OPTIMIZING: Flexible bead movers =', len(dof.get_floppy_body_movers()))
mc.optimize(10) #Note:  decreased the number of shuffle iterations
IMP.rmf.save_frame(f)
del f
# SHUFFLING -----------------------------------------------

print('SHUFFLING TIME: ', round(time.time() - start_time, 2))
if '--shuffle_only' in sys.argv:
    quit(0)

# EM Restraint
density_particles = IMP.atom.Selection(root_hierarchy, molecule='PKP1a',
                                       representation_type=IMP.atom.DENSITIES)
density_particles2 = IMP.atom.Selection(root_hierarchy, molecule='GPKP1a',
                                        representation_type=IMP.atom.DENSITIES)
density_particles = (density_particles | density_particles2).get_selected_particles()
em_kwargs = {'target_fn': './data/gmm/' + dens[0], 'slope': EM_SLOPE_PKP, 'scale_target_to_mass': True,
             'weight': EM_WEIGHT_PKP}
em_restraint = IMP.pmi.restraints.em.GaussianEMRestraint(density_particles, **em_kwargs)
lab = "PKPStructureGPKP"
em_restraint.set_label(lab)
em_restraint.add_to_model()
restraint_list.append(em_restraint)

density_particles = IMP.atom.Selection(root_hierarchy, molecules=['PG', 'DP', 'DSC1', 'DSG1'],
                                       representation_type=IMP.atom.DENSITIES)
lab = "fullPGDP"
em_kwargs = {'target_fn': './data/gmm/' + dens[1], 'slope': EM_SLOPE_PGDP, 'scale_target_to_mass': True,
             'weight': EM_WEIGHT_PGDP}
em_restraint = IMP.pmi.restraints.em.GaussianEMRestraint(density_particles.get_selected_particles(), **em_kwargs)
em_restraint.set_label(lab)
em_restraint.add_to_model()
restraint_list.append(em_restraint)

# ImmunoGold restraint
# PKP
for i in range(PKP_NUMBER):
    selection_tuple = (1, 285, 'PKP1a', i, None)
    plist = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, SAMGR_RES)
    selection_tuple2 = (286, 726, 'PKP1a', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple2, SAMGR_RES)
    yGaussian = SingleAxisMinGaussianRestraintC(mdl, plist, 115 + 158, 11, 1, 'YGaussianNTermPKP' + str(i), SAMGR_NPKP)
    yGaussian.add_to_model()
    restraint_list.append(yGaussian)
    yGaussian = SingleAxisMinGaussianRestraintC(mdl, plist2, 115 + 42, 11, 1, 'YGaussianCTermPKP' + str(i), SAMGR_CPKP)
    yGaussian.add_to_model()
    restraint_list.append(yGaussian)
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

# Binding Restraints
# PKP1 - DSG1
for i in range(len(dsg_molecules)):
    selection_tuple = (70, 213, 'PKP1a', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (570, 842, 'DSG1', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DSG1_PKP1_MPDBR",
                                                      MPDBR_WEIGHT * Y2H_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
# PKP1 - DSC1
for i in range(len(dsc_molecules)):
    selection_tuple = (1, 747, 'PKP1a', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (715, 894, 'DSC1', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DSC1_PKP1_MPDBR",
                                                      MPDBR_WEIGHT * OVERLAY_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
for i in range(PKP_NUMBER):
    selection_tuple = (70, 213, 'PKP1a', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (570, 842, 'DSG1', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (715, 894, 'DSC1', None, None)
    plist3 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2 + plist3, 0, 1, "PKP1_DSC1G1_MPDBR",
                                                      (MPDBR_WEIGHT * (OVERLAY_FACTOR + Y2H_FACTOR)) // 2)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
# PKP1 - DP
for i in range(PKP_NUMBER):
    selection_tuple = (1, 168, 'PKP1a', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 584, 'DP', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "PKP1_DP_MPDBR",
                                                      MPDBR_WEIGHT * Y2H_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
    selection_tuple = (1, 168, 'PKP1a', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 584, 'DP', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DP_PKP1_MPDBR",
                                                      MPDBR_WEIGHT * Y2H_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
# PG - DP
for i in range(PG_NUMBER):
    selection_tuple = (1, 745, 'PG', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 584, 'DP', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "PG_DP_MPDBR",
                                                      MPDBR_WEIGHT * OVERLAY_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
    selection_tuple = (1, 745, 'PG', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 584, 'DP', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DP_PG_MPDBR",
                                                      MPDBR_WEIGHT * OVERLAY_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
# DSC1-DP
for i in range(len(dsc_molecules)):
    selection_tuple = (715, 894, 'DSC1', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 176, 'DP', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DP_DSC1_MPDBR",
                                                      MPDBR_WEIGHT * OVERLAY_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)
    selection_tuple = (715, 894, 'DSC1', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 176, 'DP', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hierarchy, selection_tuple, MPDBR_RES)
    binding_res = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "DSC1_DP_MPDBR",
                                                      MPDBR_WEIGHT * OVERLAY_FACTOR)
    binding_res.add_to_model()
    restraint_list.append(binding_res)

# Localization restraints
all_particles = IMP.atom.Selection(root_hierarchy, resolution=1000,
                                   representation_type=IMP.atom.BALLS).get_selected_particles()
cylRes = CylinderLocalizationRestraintC(mdl, all_particles, 1, 135, 135, CYLINDER_RADIUS, 10,
                                        "CylinderLocalizationMain", CYL_WEIGHT)
cylRes.add_to_model()
restraint_list.append(cylRes)
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

# -----------------------------



# Finalize the protocol output
po.finalize()
s = po.system
import ihm.dumper

with open("model_init_desmosome.cif", "w") as fh:
    ihm.dumper.write(fh, [s])

# # Dataset for EM restraint
# em, = [r for r in s.restraints
#        if isinstance(r, ihm.restraint.EM3DRestraint)]
# d = em.dataset
# print("GMM file at", d.location.path)
#
last_step = s.orphan_protocols[-1].steps[-1]
last_step.num_models_end = (
    2_250_000  # 50,000 models per run and 45 independent runs (8 cores per run)
)

protocol = po.system.orphan_protocols[-1]
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis_nmodels = 24866   # nModels in cluster0
analysis.steps.append(
    ihm.analysis.ClusterStep(
        feature="RMSD",
        num_models_begin=last_step.num_models_end,
        num_models_end=analysis_nmodels,
    )
)

mg = ihm.model.ModelGroup(name="Cluster 0")
po.system.state_groups[-1][-1].append(mg)
e = ihm.model.Ensemble(
    model_group=mg,
    num_models=analysis_nmodels,
    post_process=analysis.steps[-1],
    name="Cluster 0",
    clustering_method="Density based threshold-clustering",
    clustering_feature="RMSD",
    precision="67",
)
po.system.ensembles.append(e)

Uniprot = {
    "DP.0": "P15924",
    "DP.1": "P15924",
    "DP.2": "P15924",
    "DP.3": "P15924",

    "PKP1a.0": "Q13835-2",
    "PKP1a.1": "Q13835-2",
    "PKP1a.2": "Q13835-2",
    "PKP1a.3": "Q13835-2",

    "GPKP1a.0": "Q13835-2",
    "GPKP1a.1": "Q13835-2",
    "GPKP1a.2": "Q13835-2",

    "PG.0": "P14923",
    "PG.1": "P14923",
    "PG.2": "P14923",
    "PG.3": "P14923",

    "DSG1.0": "Q02413",
    "DSG1.1": "Q02413",

    "DSC1.0": "Q08554-1",
    "DSC1.1": "Q08554-1",
}

for prot, entry in Uniprot.items():
    ref = ihm.reference.UniProtSequence.from_accession(entry)

    #TODO with DP entry
    if prot.startswith('DP'):
        ref.alignments.append(ihm.reference.Alignment(db_begin=1,db_end=584, entity_begin=1, entity_end=584))
    if prot.startswith('GPKP1a'):
        ref.alignments.append(ihm.reference.Alignment(db_begin=244,db_end=700, entity_begin=244, entity_end=700))

    if prot.startswith('DSG1'):
        ref.alignments.append(ihm.reference.Alignment(db_begin=570,db_end=842, entity_begin=570, entity_end=842))

    if prot.startswith('DSC1'):
        ref.alignments.append(ihm.reference.Alignment(db_begin=715,db_end=894,entity_begin=715, entity_end=894))

    po.asym_units[prot].entity.references.append(ref)

m = IMP.Model()
inf1 = RMF.open_rmf_file_read_only("../../../results/main_run/cluster_center_model.rmf3")
h = IMP.rmf.create_hierarchies(inf1, m)[0]
IMP.rmf.link_hierarchies(inf1, [h])
IMP.rmf.load_frame(inf1, RMF.FrameID(0))
m.update()

model = po.add_model(e.model_group)

repo = ihm.location.Repository(
    doi="10.5281/zenodo.8035862",
    root="../../../",
    top_directory="desmosome-github",
    url="https://zenodo.org/record/13363586/files/desmosome-github.zip",
)

loc_density_list = {
    "DP.0": ["LPD_DP-N", "LPD_DP-S"],
    "DP.1": ["LPD_DP-N", "LPD_DP-S"],
    "DP.2": ["LPD_DP-N", "LPD_DP-S"],
    "DP.3": ["LPD_DP-N", "LPD_DP-S"],

    "PKP1a.0": ["LPD_PKP-N", "LPD_PKP-S","LPD_PKP-C"],
    "PKP1a.1": ["LPD_PKP-N", "LPD_PKP-S","LPD_PKP-C"],
    "PKP1a.2": ["LPD_PKP-N", "LPD_PKP-S","LPD_PKP-C"],
    "PKP1a.3": ["LPD_PKP-N", "LPD_PKP-S","LPD_PKP-C"],

    "GPKP1a.0": ["LPD_GPKP"],
    "GPKP1a.1": ["LPD_GPKP"],
    "GPKP1a.2": ["LPD_GPKP"],

    "PG.0": ["LPD_PG-N", "LPD_PG-S","LPD_PG-C"],
    "PG.1": ["LPD_PG-N", "LPD_PG-S","LPD_PG-C"],
    "PG.2": ["LPD_PG-N", "LPD_PG-S","LPD_PG-C"],
    "PG.3": ["LPD_PG-N", "LPD_PG-S","LPD_PG-C"],

    "DSG1.0": ["LPD_DSG"],
    "DSG1.1": ["LPD_DSG"],

    "DSC1.0": ["LPD_DSC"],
    "DSC1.1": ["LPD_DSC"],

}

for prot, density in loc_density_list.items():
    asym = po.asym_units[prot]
    for domain_density in density:
        loc = ihm.location.OutputFileLocation(
            "../../../results/main_run/" + domain_density + ".mrc"
        )
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
        e.densities.append(den)

po.system.update_locations_in_repositories([repo])
po.finalize()

with open("model_desmosome.cif", "w") as fh:
    ihm.dumper.write(fh, [po.system])

import ihm.reader

with open("model_desmosome.cif") as fh:
    (s,) = ihm.reader.read(fh)


print('TOTAL TIME: ', round(time.time() - start_time, 2))
quit(0)
