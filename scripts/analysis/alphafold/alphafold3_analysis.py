# To get the regions of interest in a given PDB with different PAE and PLDDT

import numpy as np

import sys

from Bio.PDB  import MMCIFParser
import json
import itertools

cif = sys.argv[1] # e.g. 'dp_pkp1'
pae_json = sys.argv[2]

f = open( pae_json, "r" )
data = json.loads( f.read() )
pae = np.array( data["pae"] )


models = MMCIFParser().get_structure('cif', cif)

for model in models:

    for chain in model:
        if chain.id =="A":
            chainA = [r for r in chain.get_residues()]
        elif chain.id =="B":
            chainB = [r for r in chain.get_residues()]

    for resa,resb in itertools.product(range(len(chainA)),range(len(chainB))):

        residueA = chainA[resa]
        residueB = chainB[resb]

        if pae[resa][len(chainA)+resb] <= 5.0 or pae[len(chainA)+resb][resa] <= 5.0: # being liberal and using the cutoff on one of the PAE values for a residue pair

            if residueA["CA"].get_bfactor() > 70.0 and residueB["CA"].get_bfactor() > 70.0:

                if residueA["CA"] - residueB["CA"] < 10.0:
                    print(residueA.id[1], "A", residueB.id[1],"B")
