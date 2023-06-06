import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from scikit_posthocs import posthoc_conover, sign_plot


def check_order(main_order):  # To check that not more than 1 frame reports a 1.0 temperature at a time
    temp = []
    for i in main_order:
        temp += i.tolist()
    assert len(temp) == len(np.unique(temp)), 'Repeated 1.0 temp frame indices'
    assert max(temp) == (len(temp) - 1), 'Incorrect indices: Some indices may be missing'
    return temp


def parser(path):  # To parse all the stat files
    files = os.listdir(path)
    # Load all the appropriately named files from the folder
    stat_files = [x for x in files if re.search(r'stat[.][0-9]+[.]out', x)]
    stat_replica_files = [x for x in files if re.search(r'stat_replica[.][0-9]+[.]out', x)]
    # Check if there is no missing number in the stat files
    assert set([int(i.split('.')[1]) for i in stat_files]) == set(
        list(range(len(stat_files)))), 'Improper stat file numbering'
    for i in stat_files:  # Check if a correspondingly named replica file is present for every stat file
        x = ['stat_replica'] + i.split('.')[1:]
        assert '.'.join(x) in stat_replica_files, 'Corresponding stat_replica file missing'
    # Enough to assert the reverse of the above
    assert len(stat_files) == len(stat_replica_files), 'Corresponding stat file missing'
    # Re-creating the array to ensure that the two lists are ordered in the same order
    stat_files = ['stat.' + str(i) + '.out' for i in range(len(stat_files))]
    stat_replica_files = ['stat_replica.' + str(i) + '.out' for i in range(len(stat_files))]
    with open(path + '/' + stat_files[0]) as f:
        rd = f.readline().strip()  # The first line contains headers and environment information
        header_info = r'^.+\'STAT2HEADER_IMP_VERSIONS\'[:][ ]*([\"\'])[{].*[}]\1'
        match = re.match(header_info, rd, flags=re.DOTALL)  # Compare the expected header structure
        if not match:
            print("ERROR: Header structure different")
            quit(1)
        header_dict = eval(rd)  # Contains the number -> heading mapping as well as the environ details

    with open(path + '/' + stat_replica_files[0]) as f:
        rd = f.readline().strip()
        header_info = r'^.+\'STAT2HEADER_IMP_VERSIONS\'[:][ ]*([\"\'])[{].*[}]\1'
        match = re.match(header_info, rd, flags=re.DOTALL)
        if not match:
            print("ERROR: Header structure different")
            quit(1)
        header_dict_replica = eval(rd)  # Contains the number -> heading mapping as well as the environ details

    # build a correctly ordered trajectory
    main_array_replica = []  # contains the parsed dictionaries according to the stat_replica_files
    for i in stat_replica_files:
        with open(path + '/' + i) as f:
            rd = f.readline()  # discard line 1
            sub_array = []
            for line in f:
                sub_array.append(eval(line))
                # each element is a dictionary with integer keys
            main_array_replica.append(np.vstack([[x[y] for y in range(len(x))] for x in sub_array]))

    inverted_dict = dict()  # For heading -> number mapping
    for i in header_dict:
        if isinstance(i, int):
            inverted_dict[header_dict[i]] = i

    inverted_dict_replica = dict()  # For heading -> number mapping
    for i in header_dict_replica:
        if isinstance(i, int):
            inverted_dict_replica[header_dict_replica[i]] = i

    main_order = []  # Contains the correct order of the different frames in different replica files
    for i in main_array_replica:
        temp = i[:, inverted_dict_replica['ReplicaExchange_CurrentTemp']].flatten() == '1.0'
        main_order.append(np.where(temp)[0])

    collated_order = check_order(main_order)  # Flattened list version of main_order
    main_array = []  # contains the parsed dictionaries according to the stat_files
    for i in stat_files:
        with open(path + '/' + i) as f:
            rd = f.readline()
            sub_array = []
            for line in f:
                sub_array.append(eval(line))
            if len(sub_array) == 0:
                continue
            main_array.append(np.vstack([[x[y] for y in range(len(x))] for x in sub_array]))
    assert len(set([len(x) for x in main_array_replica])) == 1, 'Replica stat files have different number of frames'
    assert len(collated_order) == sum([len(x) for x in main_array]), 'All frames not included in main_order/array'
    assert len(collated_order) == len(main_array_replica[0]), 'All frames not included in main_order/array_replica'
    x = []
    for i in range(len(main_array)):
        x += main_array[i].tolist()
    z = sorted(zip(collated_order, x))
    z = [i[1] for i in z]  # properly ordered list of dictionaries (with stat file fields)

    x = []
    for i in range(len(main_array)):
        x += main_array_replica[i].tolist()
    z_replica = sorted(zip(collated_order, x))
    z_replica = [i[1] for i in z_replica]  # properly ordered list of dictionaries (with stat_replica file fields)
    return z, z_replica, main_array, main_array_replica, inverted_dict, inverted_dict_replica, main_order


def cutoff_convergence2(series, sigma, delta_time):
    mn = np.mean(series[-delta_time:])
    cutoff_val = mn + sigma * np.std(series[-delta_time:], ddof=1)
    return (np.array(series) < cutoff_val).tolist().index(True)


paths = ['output' + str(i + 1) for i in range(30)]

sigma = 4
delta = 2000
cutoff_foo = lambda x: cutoff_convergence2(x, sigma, delta)

ccc_array = []
total_score_array = []
evr_array = []
em_array = []
for path in paths:
    print('Plotting for ' + path)
    z, z_replica, main_array, main_array_replica, inverted_dict, inverted_dict_replica, main_order = parser(path)
    cutoff = cutoff_foo([float(i[inverted_dict['Total_Score']]) for i in z])
    conservative_cutoff = 2 * cutoff
    print(f'Cutoff: {conservative_cutoff}')

    header_list = [x for x in inverted_dict.keys() if 'Total_Score' in x]
    assert len(header_list) == 1, 'More than one Total Score variable: ' + str(header_list)
    unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
    cut_arrays = [x[conservative_cutoff:] for x in unfiltered_array_list]
    total_score_array.append(cut_arrays[0])

    header_list = [x for x in inverted_dict.keys() if 'CCC' in x]
    assert len(header_list) == 1, 'More than one CCC variable: ' + str(header_list)
    unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
    cut_arrays = [x[conservative_cutoff:] for x in unfiltered_array_list]
    ccc_array.append(cut_arrays[0])

    header_list = [x for x in inverted_dict.keys() if 'excluded' in x.lower()]
    assert len(header_list) == 1, 'More than one EVR variable: ' + str(header_list)
    unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
    cut_arrays = [x[conservative_cutoff:] for x in unfiltered_array_list]
    evr_array.append(cut_arrays[0])

    header_list = [x for x in inverted_dict.keys() if 'GaussianEMRestraint_fullPGDP' == x]
    assert len(header_list) == 1, 'More than one EM variable: ' + str(header_list)
    unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
    cut_arrays = [x[conservative_cutoff:] for x in unfiltered_array_list]
    em_array.append(cut_arrays[0])

fig, ax = plt.subplots(3, 2, figsize=(18, 18))
ax = ax.ravel()
for i in range(0, 30, 5):
    ax[i // 5].boxplot(ccc_array[i:i + 5], notch=True, sym='', bootstrap=5000)
    ax[i // 5].set_title(str(i // 5))
    ax[i // 5].set_xlabel('')
    ax[i // 5].set_ylabel('')
plt.tight_layout()
fig.savefig('plots/ccc_all_separate.png')
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
collated_array = []
collated_total_score = []
for i in range(0, 30, 5):
    collated_array.append(np.concatenate(ccc_array[i:i + 5], axis=0))
    collated_total_score.append(np.concatenate(total_score_array[i:i + 5], axis=0))
# ax.boxplot(collated_array, notch=True, sym='', bootstrap=5000, positions=[i for i in range(2, 8)], zorder=0)
ax.errorbar([i for i in range(2, 8)], [np.mean(x) for x in collated_array], yerr=[np.std(x) for x in collated_array],
            fmt='none', ecolor='black', elinewidth=3, capsize=5, capthick=3, zorder=0)
ax.scatter([i for i in range(2, 8)], [np.mean(x) for x in collated_array], s=120, color='red', marker='o', zorder=2)
ax.plot([i for i in range(2, 8)], [np.mean(x) for x in collated_array], color='black', linewidth=3, zorder=1)
ax.set_ylabel('CCC value')
fig.savefig('plots/ccc_all_combined.svg')
plt.close()

s_full, p_full = kruskal(*collated_array)
s_top10, p_top10 = kruskal(*[sorted(x)[-int(0.1 * len(x)):] for x in collated_array])
print("Full Omnibus Pval:", p_full)
print("Top10 Omnibus Pval:", p_top10)
pc_full = posthoc_conover(collated_array, p_adjust='holm')
pc_top10 = posthoc_conover([sorted(x)[-int(0.1 * len(x)):] for x in collated_array], p_adjust='holm')
print("Posthoc Pvals Full:\n", pc_full)
print("Posthoc Pvals Top10:\n", pc_top10)

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True,
                'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
sign_plot(pc_full, ax=ax, **heatmap_args)
ax.set_title('Full')
fig.savefig('plots/sign_plots_full.png')
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True,
                'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
sign_plot(pc_top10, ax=ax, **heatmap_args)
ax.set_title('Top10')
fig.savefig('plots/sign_plots_top10.png')
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter([i for i in range(2, 8)], [np.max(x) for x in collated_array], s=90, color='red', marker='+')
ax.plot([i for i in range(2, 8)], [np.max(x) for x in collated_array], color='black', linewidth=2)
fig.savefig('plots/ccc_all_combined_max.png')
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
ax.errorbar([i for i in range(2, 8)], [np.mean(sorted(x)[-int(0.1 * len(x)):]) for x in collated_array],
            yerr=[np.std(sorted(x)[-int(0.1 * len(x)):]) for x in collated_array],
            fmt='none', ecolor='black', elinewidth=3, capsize=5, capthick=3, zorder=0)
ax.scatter([i for i in range(2, 8)], [np.mean(sorted(x)[-int(0.1 * len(x)):]) for x in collated_array], s=120,
           color='red', marker='o', zorder=2)
ax.plot([i for i in range(2, 8)], [np.mean(sorted(x)[-int(0.1 * len(x)):]) for x in collated_array], color='black',
        linewidth=3, zorder=1)
fig.savefig('plots/ccc_all_combined_max_top10p.svg')
plt.close()

x = [sorted(zip(collated_total_score[i], collated_array[i])) for i in range(len(collated_array))]
plt.plot([i[0] for i in x[0]])
plt.show()
plt.close()
fig, ax = plt.subplots(figsize=(7, 7))
ms = []
ss = []
lens = []
for stoich in range(2, 8):
    temp = x[stoich - 2]
    temp = temp[:len(temp) // 10]
    lens.append(len(temp))
    only_ccc = [i[1] for i in temp]
    ms.append(np.mean(only_ccc))
    ss.append(np.std(only_ccc, ddof=1))

ax.scatter([i for i in range(2, 8)], ms, s=120, color='red', marker='o', linewidth=3, zorder=2)
ax.errorbar([i for i in range(2, 8)], ms, yerr=ss, fmt='-', elinewidth=3, linewidth=3, capthick=3, capsize=8,
            color='black', zorder=1)
# ax.errorbar([i for i in range(2, 8)], ms, yerr=[ss[i] / np.sqrt(lens[i]) for i in range(len(ss))], fmt='-',
#             capsize=5, color='black', zorder=0.5)
fig.savefig('plots/ccc_all_combined_max_top10p_score_wise.svg')
plt.close()
quit(1)
fig, ax = plt.subplots(3, 2, figsize=(13, 18))
ax = ax.ravel()
collated_array = []
evr_collated = []
em_collated = []
for i in range(0, 30, 5):
    evr = np.concatenate(evr_array[i:i + 5], axis=0)
    em = np.concatenate(em_array[i:i + 5], axis=0)
    collated_array.append(evr / em)
    evr_collated.append(evr)
    em_collated.append(em)
ax[0].boxplot(collated_array, notch=True, sym='', bootstrap=5000, positions=[i for i in range(2, 8)])
ax[0].set_title('EVR / EM')
ax[1].boxplot(evr_collated, notch=True, sym='', bootstrap=5000, positions=[i for i in range(2, 8)])
ax[1].set_title('EVR Raw')
ax[2].boxplot(em_collated, notch=True, sym='', bootstrap=5000, positions=[i for i in range(2, 8)])
ax[2].set_title('EM Raw')
ax[3].boxplot([evr_collated[i] / (i + 2) for i in range(len(evr_collated))], notch=True, sym='', bootstrap=5000,
              positions=[i for i in range(2, 8)])
ax[3].set_title('EVR Normalized (n)')
ax[4].boxplot([em_collated[i] / (i + 2) for i in range(len(em_collated))], notch=True, sym='', bootstrap=5000,
              positions=[i for i in range(2, 8)])
ax[4].set_title('EM Normalized (n)')
ax[5].boxplot([evr_collated[i] / ((i + 2) ** 2) for i in range(len(evr_collated))], notch=True, sym='',
              bootstrap=5000, positions=[i for i in range(2, 8)])
ax[5].set_title('EVR Normalized ($n^2$)')

for i in range(6):
    ax[i].set_xticks([])
fig.savefig('plots/evr_em_all_combined.png')
plt.close()
