import mrcfile
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
# compute the percentage cutoffs in the MRC range of values

inbuilt_cutoffs = [0.0167, 0.036, 0.008, 0.0117, 0.071, 0.0084, 0.0193, 0.063, 0.0035, 0.0151, 0.09]
inbuilt_names = ['DP-N', 'DP-S', 'DSC', 'DSG', 'GPKP', 'PG-C', 'PG-N', 'PG-S', 'PKP-C', 'PKP-N', 'PKP-S']

path = sys.argv[1]
all_density_files = [os.path.join(path, x) for x in os.listdir(path) if '.mrc' in x]
for d in all_density_files:
    mrc = mrcfile.open(d, 'r', permissive=True)
    n = np.sort(mrc.data.flatten())
    matching = [i for i in inbuilt_names if i in d]
    if len(matching) == 1:
        c = inbuilt_cutoffs[inbuilt_names.index(matching[0])]
    else:
        c = float(input(f'Cutoff for {d}? '))
    rng = np.linspace(np.min(n), np.max(n), 500)
    vals = [d, c, np.sum(rng < c) / 500 * 100, (np.mean(n) - c) / np.std(n), c / np.std(n)]
    print(f'{d:^20s}\t{c:.2f}\t{vals[2]:.2f}\t{vals[3]:.2f}\t{vals[4]:.2f}')
    plt.hist(n, bins=200)
    plt.yscale('log')
    y = plt.ylim()
    plt.plot([c, c], [y[0], y[1]], color='red', lw=2)
    plt.show()
