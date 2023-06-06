import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['DP-N', 'DP-S', 'DSC', 'DSG', 'GPKP', 'PG-C', 'PG-N', 'PG-S', 'PKP-C', 'PKP-N', 'PKP-S']
threshold = [0.0167, 0.036, 0.008, 0.0117, 0.071, 0.0084, 0.0193, 0.063, 0.0035, 0.0151, 0.09]
colors = ['#e31a1c', '#ef7678', '#1f78b4', '#6a3d9a', '#4daf4a', '#ffb366', '#994d00', '#ff7f00',
          '#95d293', '#377e35', '#4daf4a']

runCommand('set bgcolor white')

for i, p in enumerate(prots):
    runCommand('open LPD_' + p + '.mrc')
    runCommand('volume #' + str(i) + ' step 1 ')
    runCommand('volume #' + str(i) + ' level ' + str(threshold[i]))
    runCommand('color ' + colors[i] + ' #' + str(i))
    i += 1

# runCommand('open emd_21382.mrc')
# runCommand('volume #'+str(i)+' step 1 level 0.12 style mesh color "dark slate gray"')
