import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import tqdm


if __name__ == '__main__':
    # Contact Maps
    save_path = sys.argv[1]

    # Load coordinates
    with open('./extracted_xyzr/saved_data', 'rb') as f:
        m, p = pickle.load(f)
    with open('./contact_maps/names_indices', 'rb') as f:
        indices, residues, bead_sort, names = pickle.load(f)

    # Stahley restraint
    n_copies = 4
    axis = 1
    npkp = len(p['PKP1a']) // 4
    npg = len(p['PG']) // 4
    restraint_domains = [('PG', np.unique([indices['PG'][i] for i in range(len(indices['PG'])) if 30 <= residues['PG'][i] <= 109]))]
    restraint_domains_names = ['N-PG']
    restraint_means = [115 + 240]
    restraint_sds = [20]
    all_ds = []
    for i, name in tqdm.tqdm(enumerate(restraint_domains_names), desc='Stahley Restraints'):
        distribution_of_minimums = []
        distribution_of_averages = []
        distribution_of_signed_minimums = []
        mol, p_indices = restraint_domains[i]
        mn = restraint_means[i]
        sd = restraint_sds[i]
        n_per_copy = len(p[mol]) // n_copies
        assert len(p[mol]) % n_copies == 0, 'Wrong n_copies?'
        protein_copies = [np.array(p_indices) + i * n_per_copy for i in range(n_copies)]
        for j in tqdm.trange(m[mol].shape[0], desc=f'SAMGR {name}', smoothing=0):
            copy_domains = [m[mol][j, x, 1] - mn for x in protein_copies]
            distribution_of_signed_minimums += [x[np.argmin(np.abs(x))] for x in copy_domains]
            distribution_of_averages += [np.mean(np.abs(x)) for x in copy_domains]
            distribution_of_minimums += [np.min(np.abs(x)) for x in copy_domains]
        fig, ax = plt.subplots(2, 1, figsize=(7, 14))
        ax = ax.ravel()
        hist, bins = np.histogram(distribution_of_minimums, bins=100)
        ax[0].plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        ax[0].fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        ax[0].set_xlabel('Minimum Distance in Angstrom')
        ax[0].set_ylabel('Frequency')
        ax[0].set_title('Distribution of Minimums')
        hist, bins = np.histogram(distribution_of_averages, bins=100)
        ax[1].plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        ax[1].fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        ax[1].set_xlabel('Distance in Angstrom')
        ax[1].set_ylabel('Frequency')
        mn_mn = np.mean(distribution_of_averages)
        mn_sd = np.std(distribution_of_averages)
        ax[1].vlines([mn_mn - mn_sd, mn_mn + mn_sd], 0, ax[1].get_ylim()[1], color='red')
        ax[1].vlines([mn_mn - sd, mn_mn + sd], 0, ax[1].get_ylim()[1], color='green')
        ax[1].set_title('Distribution of Averages')
        plt.savefig(f'{save_path}/stahley_{name}_specific.png')
        plt.close(fig)
        fig = plt.figure()
        hist, bins = np.histogram(distribution_of_signed_minimums, bins=100)
        plt.plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        plt.fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        mn_mn = np.mean(distribution_of_signed_minimums)
        mn_sd = np.std(distribution_of_signed_minimums)
        plt.vlines([mn_mn - mn_sd, mn_mn + mn_sd], 0, plt.ylim()[1], color='red')
        plt.vlines([mn_mn - sd, mn_mn + sd], 0, plt.ylim()[1], color='green')
        plt.xlabel('(signed) Minimum distance in Angstrom')
        plt.ylabel('Frequency')
        plt.title('Distribution of Signed Minimums')
        plt.savefig(f'{save_path}/stahley_{name}_specific_signed.png')
        plt.close(fig)
        all_ds.append(distribution_of_signed_minimums[:])
    assert len(all_ds) == 1
    main_fig, ax_main = plt.subplots(1, 1, figsize=(10, 10))
    #ax_inset = main_fig.add_axes([0.25, 0.55, 0.3, 0.25])
    ax = [ax_main]
    vplot = ax[0].violinplot([np.array(x) for x in all_ds], showmeans=False, showextrema=False,
                             positions=np.arange(len(all_ds)))
    ax[0].plot(ax[0].get_xlim(), [0, 0], color='black')
    ax[0].plot(ax[0].get_xlim(), [- restraint_sds[0]] * 2, color='red', linestyle='-')
    ax[0].plot(ax[0].get_xlim(), [restraint_sds[0]] * 2, color='red', linestyle='-')
    restraint_domains_colors = ['#377e35', '#95d293', '#994d00', '#ffb366', '#e31a1c']
    for j, part in enumerate(vplot['bodies']):
        # part.set_facecolor(restraint_domains_colors[j])
        # part.set_edgecolor(restraint_domains_colors[j])
        part.set_facecolor('black')
        part.set_edgecolor('black')
        part.set_alpha(1)
    ax[0].set_xticks(np.arange(len(all_ds)))
    ax[0].set_xticklabels(restraint_domains_names)
    ax[0].set_ylabel('Minimum distance from the data mean Å')
    plt.savefig(f'{save_path}/stahley.svg')
    plt.close()
