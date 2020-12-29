import sys
from Bio.SeqIO.PdbIO import PdbAtomIterator

# Satwik Pasani: stochastic13 (11-Dec-2020)

# Collate two MODELLER ready .ali files into a single file to allow using multi-chain single PDB
# All the "structure" sequences from the ali files combined into one and the same for "sequence" sequences

# Input Arguments:
# 1. Set of .ali files to extract sequences from separated by comma in order
# Can enter 'PDB:' + .ali file + ':chain_id'
# The above extracts the full ATOM seq and put it as a chain in both "structure" and "sequence"
# Need to have the whole structure in the PDB (ATOM seq) covered in the alignment
# Need to have a contiguous placement of the residues in the alignment and must cover the whole of PDB
# In each ali file, one sequence and one structure must be present only
# 2. Name of the PDB file, last chain-id and the combined "sequence" protein separated by a comma
# 3. Output file path

def parse_ali_files(filename):
    if filename.split(':')[0] != 'PDB':
        with open(filename) as f:
            rd = f.read()
            rd = rd.split('>P1;')[1:3]
            rd = [x.split('\n') for x in rd]
        if rd[0][1].split(':')[0] == 'structure':
            return rd[0][2].replace('*', ''), rd[1][2].replace('*', '')
        else:
            return rd[1][2].replace('*', ''), rd[0][2].replace('*', '')
    else:
        for i in PdbAtomIterator(open(filename.split(':')[1])):
            if i.annotations['chain'] == filename.split(':')[2]:
                seq = str(i.seq).replace('X', '')
                return seq, seq


files = sys.argv[1].split(',')
seqs = [parse_ali_files(x) for x in files]
seqs_structure = '/'.join([x[0] for x in seqs])
seqs_sequence = '/'.join([x[1] for x in seqs])
pdb_name, prot_name, last_chain_id = sys.argv[2].split(',')

with open(sys.argv[3], 'w') as f:
    # output the collated alignment file (MODELLER ready)
    f.write('>P1;' + pdb_name + '\n')
    f.write('structure' + ':' + pdb_name + ':FIRST:A:LAST:' + last_chain_id + '::::\n')
    f.write(seqs_structure + '*\n')
    f.write('>P1;' + prot_name + '\n')
    f.write('sequence:::::::::\n')
    f.write(seqs_sequence + '*')
