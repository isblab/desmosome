import re
import sys
import IMP
import IMP.rmf
import IMP.atom
import RMF
import multiprocessing as mp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

import tqdm
import time
from scipy.stats import pearsonr
from scipy.spatial.distance import cdist
from collections import defaultdict
import pickle

matplotlib.use('agg')


def parse_rmf_into_radii_and_coords(location, burn_in):
    # given a rmf, extract the coordinates and radius of the particles
    m = IMP.Model()
    rmf_file = RMF.open_rmf_file_read_only(location)
    n = rmf_file.get_number_of_frames()
    print(f'Number of Models in the file {location} is {n} with the initial {burn_in} removed')
    h = IMP.rmf.create_hierarchies(rmf_file, m)[0]
    model_wise_coords = []  # all particle coordinates
    molecule_wise_coords = []
    # to be converted later to a dictionary with molecule name -> models x particle coordinates
    particle_radii = []  # all particle radii
    molecule_wise_radii = defaultdict(list)  # molecule name -> radius of all the particles in order
    molecule_wise_names = defaultdict(list)  # molecule name -> names of all the particles in order
    for i in tqdm.trange(burn_in, n):
        IMP.rmf.load_frame(rmf_file, i)
        m.update()
        all_particles = [x.get_particle() for x in IMP.atom.get_leaves(h)]
        all_coords = []
        for p in all_particles:
            xyzr = IMP.core.XYZR(p)
            all_coords.append(np.array(xyzr.get_coordinates()))
            if i == burn_in:
                particle_radii.append(xyzr.get_radius())
        molecules = h.get_children()[0].get_children()
        molecule_names = [x.get_name() for x in molecules]
        molecule_coords = defaultdict(list)
        # TODO: Use a single pass over particles to populate it molecule wise too
        for name, mol in zip(molecule_names, molecules):
            sub_array = []
            radii = []
            for x in IMP.atom.get_leaves(mol):
                xyzr = IMP.core.XYZR(x.get_particle())
                sub_array.append(np.array(xyzr.get_coordinates()))
                radii.append(xyzr.get_radius())
            molecule_coords[name] += sub_array
            if i == burn_in:
                molecule_wise_radii[name] += radii[:]
                molecule_wise_names[name] += [x.get_name() for x in IMP.atom.get_leaves(mol)]
        model_wise_coords.append(all_coords)
        molecule_wise_coords.append(molecule_coords)
    model_wise_coords = np.array(model_wise_coords)
    molecule_wise_coords_dict = dict()
    for i in molecule_wise_coords[0]:
        molecule_wise_coords_dict[i] = np.array([j[i] for j in molecule_wise_coords])
    print(f'Finished loading data from {location}')
    return model_wise_coords, molecule_wise_coords_dict, molecule_wise_names, molecule_wise_radii, particle_radii


def autocorrelation(matrix, step_sequence=None):
    # matrix can be specific to a single molecule or all the particles
    if step_sequence is None:
        step_sequence = np.arange(100)
    n_models = matrix.shape[0]
    n_particles = matrix.shape[1]
    # separately calculate the correlation across the three coordinate axes
    rhox_all = []
    rhoy_all = []
    rhoz_all = []
    px_all = []
    py_all = []
    pz_all = []
    for step in step_sequence:
        rhox_step = []
        rhoy_step = []
        rhoz_step = []
        px_step = []
        py_step = []
        pz_step = []
        for particle in range(n_particles):
            rhox, px = pearsonr(matrix[:n_models - step, particle, 0], matrix[step:, particle, 0])
            rhoy, py = pearsonr(matrix[:n_models - step, particle, 1], matrix[step:, particle, 1])
            rhoz, pz = pearsonr(matrix[:n_models - step, particle, 2], matrix[step:, particle, 2])
            rhox_step.append(rhox)
            px_step.append(px)
            rhoy_step.append(rhoy)
            py_step.append(py)
            rhoz_step.append(rhoz)
            pz_step.append(pz)
        rhox_all.append(rhox_step)
        rhoy_all.append(rhoy_step)
        rhoz_all.append(rhoz_step)
        px_all.append(px_step)
        py_all.append(py_step)
        pz_all.append(pz_step)
    return [np.array(x) for x in [rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all]]


def sort_only_columns(matrix):
    # to keep the columns intact and sort them
    return np.array(sorted(matrix.T, key=lambda x: x.tolist())).T


def plot_autocorrel_per_rmf(rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all, name, save_path, step_sequence):
    # plot all the autocorrelation decay plots for a single set of rho/p-val matrices
    if step_sequence is None:
        step_sequence = np.arange(100)
    for which_coord, rho_all, p_all in zip(['x', 'y', 'z'], [rhox_all, rhoy_all, rhoz_all], [px_all, py_all, pz_all]):
        fig, ax = plt.subplots(figsize=(10, 10))
        for i in range(rho_all.shape[1]):
            ax.plot(step_sequence, rho_all[:, i], color='black', alpha=0.1, zorder=0)
        uq = np.unique(rho_all)
        if len(uq) == 1:
            print(f'Skipping {which_coord} coordinate in {name} because it has only one value: {uq}')
        else:
            ax.plot(step_sequence, np.nanmean(rho_all, axis=1).flatten(), color='red', lw=2, zorder=100)
            ax.set_title(f'Rho_{which_coord}')
            fig.savefig(f'{save_path}/{name}_autocorrel_rho{which_coord}.png')
            if len(np.unique(p_all)) > 1:
                fig, ax = plt.subplots(figsize=(10, 10))
                non_zero_min = np.min(p_all[p_all > 0])
                p_all[p_all == 0] = non_zero_min
                xx = ax.imshow(np.log(sort_only_columns(p_all)))
                plt.colorbar(xx, ax=ax)
                ax.set_title(f'pval_{which_coord}')
                fig.savefig(f'{save_path}/{name}_autocorrel_pval{which_coord}.png')
            else:
                print(f'All {which_coord} pvals are the same value (probably zero) for {name}')
        plt.close('all')


def foo_autocorrel_worker(r, ind, save_path, step_sequence, burn_in, q):
    m1, m2, names, p1, p2 = parse_rmf_into_radii_and_coords(r, burn_in)
    rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all = autocorrelation(m1, step_sequence)
    plot_autocorrel_per_rmf(rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all, f'{ind}_all', save_path,
                            step_sequence)
    mol_dict = dict()
    for i in m2:
        a, b, c, d, e, f = autocorrelation(m2[i], step_sequence)
        plot_autocorrel_per_rmf(a, b, c, d, e, f, f'{ind}_{i}', save_path, step_sequence)
        mol_dict[f'{i}_x'] = a
        mol_dict[f'{i}_y'] = b
        mol_dict[f'{i}_z'] = c
    q.put((rhox_all, rhoy_all, rhoz_all, mol_dict, names))
    return True


def foo_rmf_parser(r, ind, burn_in, q_save):
    m1, m2, names, p1, p2 = parse_rmf_into_radii_and_coords(r, burn_in)
    q_save.put((m2, p1, ind))
    return names


def foo_saver(total, q, save_path):
    m_aggregate = None
    p_aggregate = None
    for i in range(total):
        m, p, name = q.get()
        if p_aggregate is None:
            p_aggregate = p
        else:
            for k in p:
                n = len(p[k])
                assert all([p[k][j] == p_aggregate[k][j] for j in range(n)]), f'Radii are inconsistent: {i}'
        if m_aggregate is None:
            m_aggregate = m
        else:
            for k in m_aggregate:
                m_aggregate[k] = np.vstack([m_aggregate[k], m[k]])
    with open(f'{save_path}/saved_data', 'wb') as f:
        pickle.dump([m_aggregate, p_aggregate], f)
    print('Data saved')


def foo_autocorrel_averager(q, qf, total):
    current_avg = None
    names = None
    old_total = total
    tq = tqdm.tqdm(total=old_total, desc='ac_progress', smoothing=0)
    while total > 0:
        res = q.get()
        if current_avg is None:
            current_avg = [x / old_total for x in list(res)[:3]] + [res[3]]
            names = res[-1]
            for k in res[3]:
                current_avg[3][k] /= old_total
        else:
            for i in range(3):
                current_avg[i] += res[i] / old_total
            for k in res[3]:
                current_avg[3][k] += res[3][k] / old_total
        total = total - 1
        tq.update()
    tq.close()
    qf.put((current_avg, names))


def foo_contact_maps_worker(path, molecule_pair, indices, residues, save_path, contact_cutoff=10):
    with open(f'{path}', 'rb') as f:
        m, p = pickle.load(f)
    n_models = m[list(m.keys())[0]].shape[0]
    mol1, mol2 = molecule_pair
    n_part1 = m[mol1].shape[1]
    n_part2 = m[mol2].shape[1]
    n_fragments1 = len(np.unique(indices[mol1]))
    n_fragments2 = len(np.unique(indices[mol2]))
    overall_dist_matrix = np.zeros((n_fragments1, n_fragments2), dtype=float)
    overall_dist_matrix_adjusted = np.zeros((n_fragments1, n_fragments2), dtype=float)
    overall_dist_matrix_boolean = np.zeros((n_fragments1, n_fragments2), dtype=float)
    for model in range(n_models):
        dat1 = m[mol1][model, :, :]
        dat2 = m[mol2][model, :, :]
        dist_matrix = cdist(dat1, dat2, metric='euclidean')
        radii_sum = np.array(p[mol1])[:, np.newaxis] + np.array(p[mol2])[np.newaxis, :]
        dist_matrix_adjusted = dist_matrix - radii_sum
        dist_matrix_boolean = np.array(dist_matrix_adjusted <= contact_cutoff, dtype=np.int32)
        assert (n_part1 % n_fragments1) == 0, 'Size mismatch'
        assert (n_part2 % n_fragments2) == 0, 'Size mismatch'
        for i in range(n_fragments1):
            for j in range(n_fragments2):
                overall_dist_matrix[i, j] += np.mean(dist_matrix[i::n_fragments1, j::n_fragments2])
                overall_dist_matrix_adjusted[i, j] += np.mean(dist_matrix_adjusted[i::n_fragments1, j::n_fragments2])
                temp = dist_matrix_boolean[i::n_fragments1, j::n_fragments2].flatten()
                overall_dist_matrix_boolean[i, j] += int(np.logical_or.reduce(temp))
    overall_dist_matrix /= n_models
    overall_dist_matrix_adjusted /= n_models
    overall_dist_matrix_boolean /= n_models
    m, m_adjusted, m_boolean = (overall_dist_matrix, overall_dist_matrix_adjusted, overall_dist_matrix_boolean)
    mol1, mol2 = molecule_pair
    for matrix, plot_type in zip([m, m_adjusted, m_boolean], ['unadjusted', 'adjusted', 'boolean']):
        fig, ax = plt.subplots(figsize=(10, 10))
        temp = ax.imshow(matrix[np.ix_(np.array(indices[mol1]), np.array(indices[mol2]))], aspect='auto')
        ax.set_xticks(np.arange(len(indices[mol2]))[::50])
        ax.set_yticks(np.arange(len(indices[mol1]))[::50])
        ax.set_xticklabels(residues[mol2][::50])
        ax.set_yticklabels(residues[mol1][::50])
        ax.set_title(f'Distance Matrix {plot_type}')
        ax.set_xlabel(mol2)
        ax.set_ylabel(mol1)
        plt.colorbar(temp)
        fig.savefig(f'{save_path}/dist_matrix_{mol1}_{mol2}_{plot_type}.png')
        plt.close('all')
    return True


def foo_contact_map_wrapper(args):
    return foo_contact_maps_worker(*args)


def process_names(names):
    indices = defaultdict(list)
    residues = defaultdict(list)
    reg = '[^0-9]*([0-9]+)[^0-9]+([0-9]+)[^0-9]*'
    f = open('running_temp_name_fragment_map.txt', 'w')
    for k in names:
        already_done = set()
        f.write(f'{k}\n')
        for i, name in enumerate(names[k]):
            if name in already_done:
                continue
            f.write(f'\t{name}: ')
            already_done.add(name)
            match = re.search(reg, name)
            assert match, f'Regex pattern failed to match {name}'
            n1 = int(match.group(1))
            n2 = int(match.group(2))
            if 'bead' in name:
                residues[k] += list(range(n1, n2 + 1))
                indices[k] += [i for _ in range(n1, n2 + 1)]
                f.write(f'{n1} -> {n2}\n')
            elif 'Fragment' in name:
                residues[k] += list(range(n1, n2))
                indices[k] += [i for _ in range(n1, n2)]
                f.write(f'{n1} -> {n2 - 1}\n')
            else:
                assert False, f'Neither bead nor fragment: {name} in {k}'
    f.close()
    return indices, residues


def get_all_data_cm(list_of_rmfs, npool=None, save_path='.', burn_in=0):
    print('Calculating only contact maps.')
    if npool is None:
        npool = min(len(list_of_rmfs), os.cpu_count())
    with mp.Pool(npool) as p:
        print('Setting up.')
        manager = mp.Manager()
        q3 = manager.Queue()
        cm_save = f'{save_path}/contact_maps'
        data_save = f'{save_path}/extracted_xyzr'
        if not os.path.isdir(data_save):
            os.mkdir(data_save)
        if not os.path.isdir(cm_save):
            os.mkdir(cm_save)
        args = [(list_of_rmfs[r], r, burn_in, q3) for r in range(len(list_of_rmfs))]
        p3 = mp.Process(target=foo_saver, args=(len(list_of_rmfs), q3, data_save), daemon=True)
        p3.start()
        print('Saving molecule-wise coordinates and radii.')
        temp = p.starmap(foo_rmf_parser, args)
        assert all(temp), 'Some of the foo_workers did not end properly while calculating parsing rmfs'
        names = temp[0]
        print('Finished!')
        indices, residues = process_names(names)
        del temp
        mol_pairs = set([(i, j) for i in names for j in names if i < j])
        args = [(f'{data_save}/saved_data', mol_pair, indices, residues, cm_save) for mol_pair in mol_pairs]
        total_successes = 0
        print('Waiting for processes to die')
        wait_time_start = time.time()
        while (time.time() - wait_time_start) < 60:
            p3.join(5)
            if not p3.is_alive():
                print('Processes have ended properly.')
                break
        if p3.is_alive():
            print('Processes have not ended. Exiting the function.')
            print(f'\tProcess p3: {p3.is_alive()}')
        print('Processing contact maps.')
        for result in tqdm.tqdm(p.imap_unordered(foo_contact_map_wrapper, args), desc='cm_mol_pairs', smoothing=0,
                                total=len(args)):
            total_successes += result
        check = (total_successes == len(mol_pairs))
        assert check, 'Some of the foo_workers did not end properly while calculating contact maps'


def get_all_data_autocorrel(list_of_rmfs, npool=None, save_path='.', step_sequence=None, burn_in=0):
    print('Calculating both autocorrelation and contact maps.')
    if npool is None:
        npool = min(len(list_of_rmfs), os.cpu_count())
    with mp.Pool(npool) as p:
        print('Setting up.')
        manager = mp.Manager()
        q1 = manager.Queue()
        q2 = manager.Queue()
        ac_save = f'{save_path}/autocorrel_plots'
        if not os.path.isdir(ac_save):
            os.mkdir(ac_save)
        args = [(list_of_rmfs[r], r, ac_save, step_sequence, burn_in, q1) for r in range(len(list_of_rmfs))]
        p2 = mp.Process(target=foo_autocorrel_averager, args=(q1, q2, len(list_of_rmfs)), daemon=True)
        p2.start()
        print('Processing autocorrelation data and saving molecule-wise coordinates and radii.')
        temp = p.starmap(foo_autocorrel_worker, args)
        assert all(temp), 'Some of the foo_workers did not end properly while calculating autocorrelation'
        averaged_result, names = q2.get()
        print('Finished!')
        print('Waiting for processes to die')
        wait_time_start = time.time()
        while (time.time() - wait_time_start) < 60:
            p2.join(5)
            if not p2.is_alive():
                print('Processes have ended properly.')
                break
        if p2.is_alive():
            print('Processes have not ended. Exiting the function.')
            print(f'\tProcess p2: {p2.is_alive()}')
    print('Plotting averaged plot')
    plot_autocorrel_per_rmf(averaged_result[0], averaged_result[1], averaged_result[2], [0], [0], [0], f'_averaged_all',
                            ac_save,
                            step_sequence)
    print('Plotting averaged plot molecule wise')
    for k in np.unique(['_'.join(x.split('_')[:-1]) for x in averaged_result[-1]]):
        plot_autocorrel_per_rmf(averaged_result[-1][f'{k}_x'], averaged_result[-1][f'{k}_y'],
                                averaged_result[-1][f'{k}_z'], [0],
                                [0], [0], f'_averaged_{k}', ac_save, step_sequence)
    print('Writing particle wise autocorrelations')
    with open(f'{ac_save}/_averaged_per_particle_autocorrel.csv', 'w') as f:
        headings = ','.join([f'{a}_{n}' for a in list('xyz') for n in range(95, 0, -5)])
        f.write(f'mol_name,bead_name,{headings}')
        for k in names:
            for p in range(len(names[k])):
                f.write(f'\n{k},{names[k][p]}')
                for coord in list('xyz'):
                    for perc in range(95, 0, -5):
                        a = np.array(averaged_result[-1][f'{k}_{coord}'][:, p].flatten(), dtype=float)
                        x = (a < (perc / 100)).tolist()
                        if True in x:
                            ind = x.index(True)
                        else:
                            ind = 'NA'
                        f.write(f',{ind}')


if __name__ == '__main__':
    # sys.argv -> location of all the output dirs, step sequence arguments ('none' for only cm),
    # plot saving path, number of cores, burn-in n-frames, (for cm only) cluster number
    save_path = sys.argv[3]
    n_pool = int(sys.argv[4])
    burn_in = int(sys.argv[5])
    if sys.argv[2] == 'none':
        cluster_number = int(sys.argv[6])
        location = [x for x in os.listdir(sys.argv[1]) if re.search(f'[AB]_gsm_clust{cluster_number}.rmf3', x)]
        rmfs = [f'{sys.argv[1]}/{x}' for x in location]
        get_all_data_cm(rmfs, n_pool, save_path, burn_in)
    else:
        location = [x for x in os.listdir(sys.argv[1]) if
                    ('output' == x[:6]) and os.path.isdir(f'{sys.argv[1]}/{x}') and (
                            len(os.listdir(f'{sys.argv[1]}/{x}')) > 0)]
        rmfs = []
        for x in location:
            sub_rmfs = os.listdir(f'{sys.argv[1]}/{x}/rmfs')
            sub_rmfs = [f'{sys.argv[1]}/{x}/rmfs/{i}' for i in sub_rmfs if i[-5:] == '.rmf3']
            rmfs += sub_rmfs
        step_seq = np.arange(*list(map(int, sys.argv[2].split(':'))))
        get_all_data_autocorrel(rmfs, n_pool, save_path, step_seq, burn_in)
