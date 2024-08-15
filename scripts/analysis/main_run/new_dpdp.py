import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import tqdm
from scipy.spatial.distance import cdist


if __name__ == '__main__':
    save_path = sys.argv[1]

    # Load coordinates
    with open('./extracted_xyzr/saved_data', 'rb') as f:
        m, p = pickle.load(f)
    with open('./contact_maps/names_indices', 'rb') as f:
        indices, residues, bead_sort, names = pickle.load(f)

    n_copies = {'PKP1a': 4, 'DP': 4, 'PG': 4, 'DSG1': 2, 'DSC1': 2}
    all_ds = []
    # Binding Validation data for DP - DP
    d1, len1, offset1 = ('DP', 23, 0)
    d2, len2, offset2 = ('DP', 23, 0)  # later, the actual value (DP, 9, 0) will be accounted for
    n_d1 = len(p[d1]) // n_copies[d1]
    n_d2 = len(p[d2]) // n_copies[d2]
    assert len(p[d1]) % n_copies[d1] == 0
    assert len(p[d2]) % n_copies[d2] == 0
    assert n_d1 == n_d2 == len1 == len2  # This is only to make later code easier
    protein_copies1 = [np.arange(offset1, offset1 + len1) + i * n_d1 for i in range(n_copies[d1])]
    protein_copies2 = [np.arange(offset2, offset2 + len2) + i * n_d2 for i in range(n_copies[d2])]
    protein_copy_coords1 = [m[d1][:, x, :] for x in protein_copies1]
    protein_copy_coords2 = [m[d2][:, x, :] for x in protein_copies2]
    protein_bead_sizes1 = np.concatenate([np.array(p[d1])[x] for x in protein_copies1])
    protein_bead_sizes2 = np.concatenate([np.array(p[d2])[x] for x in protein_copies2])
    distribution_of_minimums = []
    # distribution_of_averages = []  # column averages or row averages or both?
    for j in tqdm.trange(m[d1].shape[0], desc='DP-DP', smoothing=0):
        c1 = [x[j, :, :] for x in protein_copy_coords1]
        c2 = [x[j, :, :] for x in protein_copy_coords2]
        dist_matrix = cdist(np.concatenate(c1, axis=0), np.concatenate(c2, axis=0), metric='euclidean')
        # adjust for the self-copy inclusion
        self_distance_mask = dist_matrix < 1e-5  # zero distance is not expected elsewhere
        assert np.sum(self_distance_mask.flatten()) == (n_d1 * n_copies['DP'])
        for i in range(n_copies['DP']):
            region_of_overlap = np.arange(i * n_d1, i * n_d1 + n_d1)
            for ind1 in region_of_overlap:
                for ind2 in region_of_overlap:
                    dist_matrix[ind1, ind2] = np.inf
        radii_sum = protein_bead_sizes1[:, np.newaxis] + protein_bead_sizes2[np.newaxis, :]
        dist_matrix_adjusted = dist_matrix - radii_sum
        dist_matrix_adjusted[dist_matrix_adjusted < 0] = 0  # to consider overlapping beads as contacts
        # Now, take only the part of the distmatrix for (DP, 9, 0)
        indices_to_pick = np.array(list(range(9)))
        indices_to_pick = [indices_to_pick + i * n_d1 for i in range(n_copies['DP'])]
        indices_to_pick = np.concatenate(indices_to_pick)
        dist_matrix_adjusted = dist_matrix_adjusted[:, indices_to_pick]
        distribution_of_minimums.append(np.min(dist_matrix_adjusted.flatten()))
        # distribution_of_averages.append(np.mean(dist_matrix_adjusted.flatten()))
    all_ds.append(distribution_of_minimums)
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    # remove outliers only if all values are not zero (or all values are not the same)
    all_ds = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] if len(np.unique(x)) > 5 else x for x in all_ds]
    vplot = ax.violinplot(all_ds, positions=np.arange(len(all_ds)), showmeans=False, showextrema=False)
    for j, part in enumerate(vplot['bodies']):
        part.set_color('black')
        part.set_alpha(1)
        for perc in np.percentile(np.array(all_ds[j]), [5, 50, 95]):
            ax.plot([j - 0.25, j + 0.25], [perc, perc], color='red', zorder=100)
    ax.set_xticks(np.arange(len(all_ds)))
    ax.set_xticklabels(['DP-DP'])
    ax.set_xlabel('Binding Validation Data')
    ax.set_ylabel('Distribution of model-wise minima')
    plt.savefig(f'{save_path}/binding_all_validation.svg')
    plt.close()
