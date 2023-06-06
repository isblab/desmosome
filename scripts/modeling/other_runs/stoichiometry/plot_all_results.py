import os
import re
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt


def check_order(main_order):  # To check that not more than 1 frame reports a 1.0 temperature at a time
    temp = []
    for i in main_order:
        temp += i.tolist()
    assert len(temp) == len(np.unique(temp)), 'Repeated 1.0 temp frame indices'
    assert max(temp) == (len(temp) - 1), 'Incorrect indices: Some indices may be missing'
    return temp


def sort_the_replica_exchanges_lowest_temp(main_order):
    m = []  # To flatten the main_order list
    replica_index = []  # Integer indices given to each replica
    for i in range(len(main_order)):
        m += main_order[i].tolist()
        replica_index += ([i for x in range(len(main_order[i]))])
    sorted_replicas = sorted(zip(m, replica_index))  # Sort the indices according the main_order
    sorted_replicas = [x[1] for x in sorted_replicas]  # extract the indices
    x = [False] + (np.diff(sorted_replicas) != 0).tolist()  # Boolean array marking exchanges
    return sorted_replicas, x


def sort_the_replica_exchanges_all_temp(main_array_replica, inverted_dict_replica):
    x = np.array([0 for i in range(len(main_array_replica[0]))], dtype=np.int32)  # Total exchanges at each step
    for i in main_array_replica:
        # Get the temperatures for each frame for this particular replica
        curr_replica_temp = i[:, inverted_dict_replica['ReplicaExchange_CurrentTemp']].flatten()
        curr_replica_temp = [float(x) for x in curr_replica_temp]
        # Check whenever the temperature changes indicating an exchange
        curr_replica_exchanges = [0] + (np.abs(np.diff(curr_replica_temp)) > 1e-5).tolist()
        curr_replica_exchanges = np.array(curr_replica_exchanges, dtype=np.int32)
        x += curr_replica_exchanges
    return x / 2  # Since each exchange is added twice in the two exchanging replicas


def plot_exchanges(ax, exchanges, all_exchanges, rectangle_size=0.03, width=3):
    ylim = ax.get_ylim()  # To put the exchanges heat plot
    padding = 0.25 * rectangle_size * (ylim[1] - ylim[0])  # To pad between the two rectangles
    rectangle_size = rectangle_size * (ylim[1] - ylim[0])  # Percentage of the total width
    n_points = int(6 * width)  # 3 sd on either side
    # Compute a gaussian convolution with the given SD (= width) over the number of exchanges
    gaussian_kernel = np.exp(-0.5 * ((np.linspace(-n_points // 2, n_points // 2, n_points) / width) ** 2))
    # plot the minimum temperature exchanges first
    exchanges = np.array(exchanges, dtype=np.float64)
    values = np.convolve(exchanges, gaussian_kernel, mode='same')
    values = np.array([values.tolist(), values.tolist()])  # For ease of plotting below
    ax.pcolormesh([i for i in range(len(exchanges))], [ylim[1] + rectangle_size, ylim[1]], values, shading='gouraud',
                  cmap='magma')
    # plot all the exchanges
    exchanges = np.array(all_exchanges, dtype=np.float64)
    values = np.convolve(exchanges, gaussian_kernel, mode='same')
    values = np.array([values.tolist(), values.tolist()])  # For ease of plotting below
    ax.pcolormesh([i for i in range(len(exchanges))],
                  [ylim[1] + 2 * rectangle_size + padding, ylim[1] + rectangle_size + padding], values,
                  shading='gouraud', cmap='magma')


def correct_mc_cumulative(mc_array, min_temp_exchanges, num_steps=10):
    # Cumulative number of MonteCarlo steps
    adjusted_mc_array = []
    for sub_mc in mc_array:
        number_of_steps = np.arange(1, len(sub_mc) + 1) * num_steps
        corrected_array = [sub_mc[0] * num_steps]
        # contains the "actual" number of accepted steps per inter-frame interval
        for i in range(1, len(sub_mc)):
            nsteps = sub_mc[i] * number_of_steps[i] - sub_mc[i - 1] * number_of_steps[i - 1]
            corrected_array.append(nsteps)
        # At all exchanges, the denominator changes, and hence, the numbers are not valid
        corrected_array = np.array(corrected_array) / num_steps
        corrected_array[np.array(min_temp_exchanges)] = np.NAN
        adjusted_mc_array.append(corrected_array)
    return adjusted_mc_array


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

mc_arrays = defaultdict(list)
labs = ['dp_structure', 'pg_dsc1_structure', 'dp_srb', 'pg_dsc1_srb', 'BallMover']
labs2 = ['dp', 'pg', 'dp_srb', 'pg_srb', 'bead']
for path in paths:
    print('Plotting for ' + path)
    z, z_replica, main_array, main_array_replica, inverted_dict, inverted_dict_replica, main_order = parser(path)
    sorted_replicas, exchange_indices = sort_the_replica_exchanges_lowest_temp(main_order)
    all_exchange_indices = sort_the_replica_exchanges_all_temp(main_array_replica, inverted_dict_replica)
    cutoff = cutoff_foo([float(i[inverted_dict['Total_Score']]) for i in z])
    conservative_cutoff = 2 * cutoff
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot([float(i[inverted_dict['Total_Score']]) for i in z][cutoff:], linestyle=':', color='black', linewidth=3,
            label='Min Temp')
    index = 0
    for i in main_array_replica:
        ax.plot(np.array(i[:, inverted_dict_replica['score']].flatten(), dtype=float)[cutoff:], alpha=0.5,
                label=str(index))
        index += 1
    ax.set_xlabel('Frame Numbers after cutoff: ' + str(cutoff))
    ax.set_ylabel('Total Score')
    ax.set_title('Total Score over time')
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='x-small')
    plt.subplots_adjust(right=0.8, left=0.1)
    # Plot the replica exchanges
    plot_exchanges(ax, exchange_indices[cutoff:], all_exchange_indices[cutoff:],
                   width=len(exchange_indices[cutoff:]) // 1000)
    fig.savefig('plots/total_score_cutoff_' + path + '.png')
    plt.close()
    # unfiltered plot only for a rough estimate of the filtering amount
    fig, ax = plt.subplots()
    ax.plot([float(i[inverted_dict['Total_Score']]) for i in z], linestyle=':', color='black', linewidth=3)
    for i in main_array_replica:
        ax.plot(np.array(i[:, inverted_dict_replica['score']].flatten(), dtype=float), alpha=0.5)
    ax.set_xlabel('Frame Numbers')
    ax.set_ylabel('Total Score')
    ax.set_title('Total Score over time')
    ax.axvline(cutoff, color='black')
    ax.fill_between([0, cutoff], [0, 0], [ax.get_ylim()[1]] * 2, color='red', alpha=0.1)
    fig.savefig('plots/total_score_raw_' + path + '.png')
    plt.close()

    for a, b in zip(labs, labs2):
        header_list = [x for x in inverted_dict.keys() if ('MonteCarlo_Acceptance_' in str(x)) and (a in str(x))]
        unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
        adjusted_unfiltered_array_list = correct_mc_cumulative(unfiltered_array_list, exchange_indices)
        cut_arrays = [x[conservative_cutoff:] for x in adjusted_unfiltered_array_list]
        flattened_array = np.mean(np.array(cut_arrays, dtype=float), axis=0)
        mc_arrays[b].append(np.nanmean(flattened_array))

    cut_indices = all_exchange_indices[conservative_cutoff:]
    global_average = np.sum(cut_indices) / len(cut_indices)
    cut_indices2 = exchange_indices[conservative_cutoff:]
    global_average2 = np.sum(cut_indices2) / len(cut_indices2)
    mc_arrays['rex_all_temp'].append(global_average / 4)
    mc_arrays['rex_min_temp'].append(global_average2)

for b in labs2:
    fig, ax = plt.subplots(figsize=(7, 7))
    vec = mc_arrays[b]
    vec_means = [np.nanmean(vec[i:i + 5]) for i in range(0, 30, 5)]
    vec_std = [np.nanstd(vec[i:i + 5], ddof=1) for i in range(0, 30, 5)]
    ax.errorbar([i for i in range(2, 8)], vec_means, yerr=vec_std, fmt='r.', ecolor='black', elinewidth=2, capsize=5,
                capthick=2)
    ax.set_xlabel('Stoichiometry')
    ax.plot([1, 8], [0.3, 0.3], color='red', linestyle='--')
    ax.plot([1, 8], [0.6, 0.6], color='red', linestyle='--')
    ax.set_ylabel('MC Acceptance Ratios for ' + b)
    fig.savefig('plots/mc_' + b + '_all.png')
    plt.close()

for b in ['rex_all_temp', 'rex_min_temp']:
    fig, ax = plt.subplots(figsize=(7, 7))
    vec = mc_arrays[b]
    vec_means = [np.mean(vec[i:i + 5]) for i in range(0, 30, 5)]
    vec_std = [np.std(vec[i:i + 5], ddof=1) for i in range(0, 30, 5)]
    ax.errorbar([i for i in range(2, 8)], vec_means, yerr=vec_std, fmt='r.', ecolor='black', elinewidth=2, capsize=5,
                capthick=2)
    ax.set_xlabel('Stoichiometry')
    ax.plot([1, 8], [0.15, 0.15], color='red', linestyle='--')
    ax.plot([1, 8], [0.25, 0.25], color='red', linestyle='--')
    ax.set_ylabel('REX exchange ratios for ' + b)
    fig.savefig('plots/rex_' + b + '_all.png')
    plt.close()
