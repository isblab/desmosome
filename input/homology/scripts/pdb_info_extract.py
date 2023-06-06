import sys
import os
import re
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqIO.PdbIO import PdbSeqresIterator, PdbAtomIterator

# Satwik Pasani: stochastic13 (14-Dec-2020)

# Extracts residue numbers

# Input Arguments (separated by space):

# 1. PDB path of the chosen template
# 2. Path to Output directory (don't use current dir; no trailing slashes; use forward slash '/')

# prot_index: index of the residue in the protein
# seq_index: index of the residue in a sequence (the full SEQRES sequence unless mentioned)
prot_index = {}
seq_index = {}
captured_hetatm = ['SEP', 'MSE']  # the atoms for which the MODELLER ali file will have a '.' residue
structure = PDBParser(QUIET=True).get_structure('main', sys.argv[1])
for chain in structure[0]:
    prot_index['structure_begin'] = [x for x in structure[0][chain.id].get_residues()][0].id[1]
    seqres = PdbSeqresIterator(open(sys.argv[1]))  # Reads from the SEQRES records (SEP -> S)
    atomres = PdbAtomIterator(open(sys.argv[1]))  # Reads from the ATOM records, (SEP -> S; missing -> X)
    for i in seqres:
        if i.annotations['chain'] == chain.id:  # Choose the chain
            full_seq = str(i.seq)
    for i in atomres:
        if i.annotations['chain'] == chain.id:  # Choose the chain
            structure_seq = str(i.seq)
    id_start = re.finditer(structure_seq.replace('X', '.'), full_seq)  # regex search of the ATOM seq on SEQRES seq
    count = 0
    for x in id_start:
        seq_index['structure_begin'] = x.start(0)
        end = x.end(0)
        count += 1
    assert count == 1, 'More than one match found/No matches found.\n' + full_seq + '\n' + structure_seq

    if os.path.isfile(sys.argv[2] + '/pdb_res_num_table.txt'):
        f = open(sys.argv[2] + '/pdb_res_num_table.txt', 'a')
    else:
        f = open(sys.argv[2] + '/pdb_res_num_table.txt', 'w')
        f.write('pdb\tprot_index_SEQRES\tSEQRES_len\tprot_index_ATOM\tATOM_num\tATOM_len\tATOM_missing\n')
    f.write(sys.argv[1] + '_' + chain.id + '\t')
    f.write(str(prot_index['structure_begin'] - seq_index['structure_begin']) + '\t')
    f.write(str(len(full_seq)) + '\t')
    f.write(str(prot_index['structure_begin']) + '\t')
    f.write(str(len(structure_seq) - list(structure_seq).count('X')) + '\t')
    f.write(str(len(structure_seq)) + '\t')
    all_indices = np.array(
        [x for x in range(prot_index['structure_begin'], prot_index['structure_begin'] + len(structure_seq))])
    missing_indices = all_indices[np.array(list(structure_seq)) == 'X'].tolist()
    missing_indices = [str(x) for x in missing_indices]
    f.write(';'.join(missing_indices) + '\n')
    f.close()
