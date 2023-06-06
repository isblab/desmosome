import matplotlib.pyplot as plt
import sys
from modeller import *

# Satwik Pasani: stochastic13 (23-Dec-2020)

# Plots DOPE profiles (normalized_dope and dope_hr) to identify poorly performing regions

# Input Arguments (separated by space):
# 1. Name of the PDB to analyze
# 2. '-' separated groups of residue numbers of interest to highlight, each group is ':' separated
# eg: 10:30_230:300, number x indicates the xth residue with structure


# boilerplate initialization code taken from MODELLER tutorials------------------
env = environ()
env.io.hetatm = True
env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters
# -------------------------------------------------------------------------------

mdl = model(env, file=sys.argv[1])
smdl = selection(mdl)

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax = ax.flatten()
ax[0].plot([x.energy for x in mdl.get_normalized_dope_profile()], color='red', alpha=0.6)
ax[0].plot([x.energy for x in mdl.get_normalized_dope_profile().get_smoothed(10)], color='black')
ax[1].plot([x.energy for x in smdl.get_dopehr_profile()], color='red', alpha=0.3)
ax[1].plot([x.energy for x in smdl.get_dopehr_profile().get_smoothed(10)], color='black')
regions_of_interest = [(int(x.split(':')[0]), int(x.split(':')[1])) for x in sys.argv[2].split('_')]
for i, j in regions_of_interest:
    ax[0].fill_between([i, j], [ax[0].get_ylim()[0], ax[0].get_ylim()[0]],
                       [ax[0].get_ylim()[1], ax[0].get_ylim()[1]], color='green', alpha=0.2)
    ax[1].fill_between([i, j], [ax[1].get_ylim()[0], ax[1].get_ylim()[0]],
                       [ax[1].get_ylim()[1], ax[1].get_ylim()[1]], color='green', alpha=0.2)
plt.show()
