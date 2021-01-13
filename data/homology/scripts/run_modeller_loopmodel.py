from modeller import *
from modeller.automodel import *
import sys
import os
import shutil

# Satwik Pasani: stochastic13 (14-Dec-2020)

# Run the automodel() modeling runs with enhanced MD refinement and loop modeling and a higher resolution assessment

# Input Arguments:
# 1. .ali file (MODELLER ready)
# 2. "knowns", i.e. the name of the template pdb
# 3. "sequence", i.e. name of the modeled sequence
# 4. Output directory (the outputted file will be moved to the directory after the run)
# 5. PDB directory
env = environ(rand_seed=-8123)  # just stating the default seed explicitly
env.io.hetatm = True
env.io.atom_files_directory = [sys.argv[5]]
aln = alignment(env, file=sys.argv[1], align_codes=(sys.argv[2], sys.argv[3]))  # to feed in the overriden method below


class mymodel(dopehr_loopmodel):
    # overriding this method to overcome the non-selected large loop in the model (basically a debug run)
    def select_loop_atoms(self):
        return selection(self.loops(aln=aln, minlength=7, maxlength=200, insertion_ext=1, deletion_ext=1))


prior_file_list = os.listdir('.')
a = mymodel(env, alnfile=sys.argv[1], knowns=sys.argv[2], sequence=sys.argv[3],
            assess_methods=(assess.normalized_dope, assess.GA341, assess.DOPE, assess.DOPEHR),
            loop_assess_methods=(assess.normalized_dope, assess.GA341, assess.DOPE, assess.DOPEHR))
a.starting_model = 1
a.ending_model = 5
a.md_level = refine.slow

a.loop.starting_model = 1
a.loop.ending_model = 3
a.loop.md_level = refine.slow

a.make()
winner = ''
curr = 100
with open('model_summary.txt', 'w') as f:
    f.write('model_name\tmolpdf\tDOPE\tDOPEHR\tnorm_DOPE\tGA341\n')
    for i in a.outputs:
        f.write(i['name'] + '\t')
        f.write(str(i['molpdf']) + '\t')
        f.write(str(i['DOPE score']) + '\t')
        f.write(str(i['DOPE-HR score']) + '\t')
        f.write(str(i['Normalized DOPE score']) + '\t')
        f.write(str(max(i['GA341 score'])) + '\t\n')
        if float(i['Normalized DOPE score']) < curr:
            curr = float(i['Normalized DOPE score'])
            winner = i['name']
    for i in a.loop.outputs:
        f.write(i['name'] + '\t')
        f.write(str(i['molpdf']) + '\t')
        f.write(str(i['DOPE score']) + '\t')
        f.write(str(i['DOPE-HR score']) + '\t')
        f.write(str(i['Normalized DOPE score']) + '\t')
        f.write(str(max(i['GA341 score'])) + '\t\n')
        if float(i['Normalized DOPE score']) < curr:
            curr = float(i['Normalized DOPE score'])
            winner = i['name']
print('\nThe final winner: ' + winner)

for file in os.listdir('.'):
    if file not in prior_file_list:
        shutil.move(file, sys.argv[4] + '/' + file)
