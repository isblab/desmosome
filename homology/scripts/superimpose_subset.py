from Bio.PDB import PDBParser, PDBIO, Superimposer
import sys

# Satwik Pasani: stochastic13 (11-Dec-2020)

# Superimpose a PDB strcuture on another using a subset of either structures

# Input Arguments:
# 1. PDB template (i.e. fixed atoms)
# 2. PDB target (i.e. moving atoms)
# 3. start,end,chain of the residues in PDB_template (in prot_index, both inclusive)
# 4. start,end,chain of the residues in PDB_template (in prot_index, both inclusive)
# 5. Output PDB file path
# 6. Flag ("ca") to use only C-alpha atoms for alignment or use all atoms ('all')
# with 'all' flag active, different missing residues in the intervening regions will cause an error

structure_template = PDBParser(QUIET=True).get_structure('template', sys.argv[1])[0]
structure_target = PDBParser(QUIET=True).get_structure('target', sys.argv[2])[0]

selected_atoms = [[], []]
actually_appended_indices = [[], []]
start1, end1, chain1 = [x for x in sys.argv[3].split(',')]
start2, end2, chain2 = [x for x in sys.argv[4].split(',')]
start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)
assert (end1 - start1) == (end2 - start2), 'Non-matching Ranges'
use_all = None
if sys.argv[6].lower() == 'all':
    use_all = True
elif sys.argv[6].lower() == 'ca':
    use_all = False
else:
    assert False, 'Incorrect Flag'
resnames_debug = [[], []]
for res in structure_template[chain1]:
    if (int(res.id[1]) >= start1) and (int(res.id[1]) <= end1):
        resnames_debug[0].append(res.get_resname())
        actually_appended_indices[0].append(int(res.id[1]) - start1)
        if use_all:
            for atom in res.get_atoms():
                selected_atoms[0].append(atom)
        else:
            selected_atoms[0].append(res['CA'])

for res in structure_target[chain2]:
    if (int(res.id[1]) >= start1) and (int(res.id[1]) <= end1):
        resnames_debug[1].append(res.get_resname())
        actually_appended_indices[1].append(int(res.id[1]) - start2)
        if use_all:
            for atom in res.get_atoms():
                selected_atoms[1].append(atom)
        else:
            selected_atoms[1].append(res['CA'])
common_indices = [i for i in range(end1 - start1 + 1) if
                  ((i in actually_appended_indices[0]) and (i in actually_appended_indices[1]))]
selected_atoms[0] = [selected_atoms[0][i] for i in range(len(actually_appended_indices[0])) if
                     actually_appended_indices[0][i] in common_indices]
selected_atoms[1] = [selected_atoms[1][i] for i in range(len(actually_appended_indices[1])) if
                     actually_appended_indices[1][i] in common_indices]
resnames_debug[0] = [resnames_debug[0][i] for i in range(len(actually_appended_indices[0])) if
                     actually_appended_indices[0][i] in common_indices]
resnames_debug[1] = [resnames_debug[1][i] for i in range(len(actually_appended_indices[1])) if
                     actually_appended_indices[1][i] in common_indices]
assert resnames_debug[0] == resnames_debug[1], 'The Amino Acid sequences do not match'
imposer_object = Superimposer()
imposer_object.set_atoms(selected_atoms[0], selected_atoms[1])
imposer_object.apply(structure_target.get_atoms())
print('RMSD: ', imposer_object.rms)

io_object = PDBIO()
io_object.set_structure(structure_target)
io_object.save(sys.argv[5])
