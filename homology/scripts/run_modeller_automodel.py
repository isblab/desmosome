from modeller import *
from modeller.automodel import *
import sys
import os
import shutil

# Satwik Pasani: stochastic13 (11-Dec-2020)

# Run the automodel() modeling runs for a given number of models with HETATM=True

# Input Arguments:
# 1. .ali file (MODELLER ready)
# 2. "knowns", i.e. the name of the template pdb
# 3. "sequence", i.e. name of the modeled sequence
# 4. Output directory (the outputted file will be moved to the directory after the run)
# 5. PDB directory

prior_file_list = os.listdir('.')
env = environ(rand_seed=-8123)  # just stating the default seed explicitly
env.io.hetatm = True
env.io.atom_files_directory = [sys.argv[5]]
a = automodel(env, alnfile=sys.argv[1], knowns=sys.argv[2], sequence=sys.argv[3],
              assess_methods=(assess.normalized_dope, assess.GA341, assess.DOPE))
a.starting_model = 1
a.ending_model = 10
a.make()
winner = ''
curr = 100
with open('model_summary.txt', 'w') as f:
    f.write('model_name\tmolpdf\tDOPE\tnorm_DOPE\tGA341\n')
    for i in a.outputs:
        f.write(i['name'] + '\t')
        f.write(str(i['molpdf']) + '\t')
        f.write(str(i['DOPE score']) + '\t')
        f.write(str(i['Normalized DOPE score']) + '\t')
        f.write(str(max(i['GA341 score'])) + '\t\n')
        if float(i['Normalized DOPE score']) < curr:
            curr = float(i['Normalized DOPE score'])
            winner = i['name']
print('\nThe final winner: ' + winner)

for file in os.listdir('.'):
    if file not in prior_file_list:
        shutil.move(file, sys.argv[4] + '/' + file)
