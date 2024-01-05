import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.pmi.analysis
import sys
import re
import tqdm
import numpy as np
import matplotlib.pyplot as plt
    

model_file = sys.argv[1]
mdl = IMP.Model()
all_models = RMF.open_rmf_file_read_only(model_file)
hier = IMP.rmf.create_hierarchies(all_models, mdl)[0]
minimum_distance = []
for i in tqdm.trange(all_models.get_number_of_frames()):
    IMP.rmf.load_frame(all_models, RMF.FrameID(i))
    mdl.update()
    sel_dsc = IMP.atom.Selection(hier, molecule='DSC1', residue_indexes=[715]).get_selected_particles()
    sel_dsg = IMP.atom.Selection(hier, molecule='DSG1', residue_indexes=[570]).get_selected_particles()
    assert len(sel_dsc) == len(sel_dsg) == 2
    dsc_coords = [IMP.core.XYZ(p).get_coordinates() for p in sel_dsc]
    dsg_coords = [IMP.core.XYZ(p).get_coordinates() for p in sel_dsg]
    all_dists = []
    for i in dsc_coords:
        for j in dsg_coords:
            all_dists.append(np.sqrt(((np.array(i) - np.array(j)) ** 2).sum()))
    minimum_distance.append(min(all_dists))

fig, ax = plt.subplots(figsize=(10, 10))
ax.hist(minimum_distance, bins=100, color='black')
ax.set_title('Minimum DSG1-DSC1 distance for each model')
ax.set_xlabel('Distance (Å)')
ax.set_ylabel('Number of Models')
ax.set_xticks([30, 50, 70, 90, 110], [30, 50, 70, 90, 110])
ax.axvline(70, color='green', label='70 Å')
ax.axvline(np.mean(minimum_distance), color='red', label='mean')
print(f'Mean: {np.mean(minimum_distance)}\nSD: {np.std(minimum_distance)}')
ax.legend()
plt.savefig('sikora_dsc_dsg.svg')
plt.show()