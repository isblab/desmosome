import numpy as np
import os
import sys

# calculates the percentage values above the thresholds for contact maps
location = sys.argv[1]
output_location = sys.argv[2]
files = [i for i in os.listdir(location) if ('flattened_dist_matrices' in i) and ('.txt' in i)]
with open(os.path.join(output_location, 'thresholds_all_mol_pairs.txt'), 'w') as f:
    temp = 'mol_pair'
    temp2 = '\t'.join(list(map(str, np.around(np.linspace(15, 35, 21), 0))))
    f.write(f'{temp:^10}\t{temp2}\n')
for i in files:
    mol_pair = '_'.join(i.split('.')[0].split('_')[3:])
    with open(os.path.join(location, i)) as f:
        rd = f.read().strip()
        rd = rd.split('\n')
        rd = [list(map(float, i.split(','))) for i in rd if i]
        rd = np.array(rd)
    vals = []
    for c in np.linspace(0.15, 0.35, 21):
        vals.append(np.sum(rd[-1] >= c) / len(rd[-1]) * 100)
    vals = '\t'.join(list(map(str, np.around(vals, 0))))
    with open(os.path.join(output_location, 'thresholds_all_mol_pairs.txt'), 'a') as f:
        f.write(f'{mol_pair:^10}:\t{vals}\n')
