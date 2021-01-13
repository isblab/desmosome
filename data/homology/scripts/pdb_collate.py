from Bio.PDB import PDBParser, PDBIO
import os
import sys

# Satwik Pasani: stochastic13 (12-Dec-2020)

# Combine different chains from different PDBs into a single PDB

# Input Arguments:
# This combined the different chains with arbitrary numbering and may even result in an improper PDB (owing partly to
# Biopython's and the PDB file format's nuances) -Auxiliary information (eg: SEQRES) may be lost/collated improperly
# Space separated elements
# Each element is comma separated 'pdb_name','chain'
# All the entered chains are combined into a single PDB

names = [i.split(',') for i in sys.argv[1:]]
chains = []
for i, j in names[1:]:  # first name used as the final structure
    chains.append(PDBParser(QUIET=True).get_structure('temp', i)[0][j])
# Create a structure and rid it of all unnecessary chains
s = PDBParser(QUIET=True).get_structure('temp_final', names[0][0])[0]
output_pdb = PDBIO()
output_pdb.set_structure(s[names[0][1]])
output_pdb.save('temp__.pdb')

# using this round-about method instead of serial detach_child detach_parent calls because in the latter, in some PDB,
# with multiple equivalent chains, some issues crop up
s = PDBParser(QUIET=True).get_structure('final', 'temp__.pdb')[0]
for chain in s:
    chain.id = 'A'  # automatically asserts that only one chain is present

curr = 'B'
for chain in chains:
    chain.detach_parent()
    chain.id = curr
    s.add(chain)
    curr = chr(ord(curr) + 1)

output_pdb = PDBIO()
output_pdb.set_structure(s)
output_pdb.save('combined.pdb')
os.remove('temp__.pdb')
