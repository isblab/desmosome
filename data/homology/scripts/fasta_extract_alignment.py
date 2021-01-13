import sys
import os
import re
from Bio.PDB import PDBParser
from Bio.SeqIO.PdbIO import PdbSeqresIterator, PdbAtomIterator
from Bio.PDB.Polypeptide import three_to_one

# Satwik Pasani: stochastic13 (23-Dec-2020)

# Creates MODELLER ready alignment files from the SEQRES sequence of PDB files 
# (in cases where there are missing residues in a PDB and we need to model them)

# Input Arguments (separated by space):
# 1. PDB path of the chosen template
# 2. Chain of the PDB to use
# 3. Path to Output directory (don't use current dir; no trailing slashes; use forward slash '/')
# 4. Name of the query to store in the alignment

# prot_index: index of the residue in the protein
# seq_index: index of the residue in a sequence (the full SEQRES sequence unless mentioned)

# To maintain code-uniformity, we just create a hypothetical alignment (line 37) and then use repeat-code
# from hhr_extract_alignment for downstream manipulation
prot_index = {}
seq_index = {}
captured_hetatm = ['SEP', 'MSE']  # the atoms for which the MODELLER ali file will have a '.' residue
structure = PDBParser(QUIET=True).get_structure('main', sys.argv[1])
prot_index['structure_begin'] = [x for x in structure[0][sys.argv[2]].get_residues()][0].id[1]
seqres = PdbSeqresIterator(open(sys.argv[1]))  # Reads from the SEQRES records (SEP -> S)
atomres = PdbAtomIterator(open(sys.argv[1]))  # Reads from the ATOM records, (SEP -> S; missing -> X)
for i in seqres:
    if i.annotations['chain'] == sys.argv[2]:  # Choose the chain
        full_seq = str(i.seq)
for i in atomres:
    if i.annotations['chain'] == sys.argv[2]:  # Choose the chain
        structure_seq = str(i.seq)
template_ali, query_ali = full_seq, full_seq
seq_index['pseudo_template_SEQRES'] = 1
seq_index['pseudo_query_input'] = 1  # assuming we queried only with the SEQRES sequence
id_start = re.finditer(structure_seq.replace('X', '.'), full_seq)  # regex search of the ATOM seq on SEQRES seq
count = 0
for x in id_start:
    seq_index['structure_begin'] = x.start(0)
    end = x.end(0)
    count += 1
assert count == 1, 'More than one match found/No matches found.\n' + full_seq + '\n' + structure_seq
final_alignment = [[], []]
iter_range = max(0, seq_index['structure_begin'] - (seq_index['pseudo_template_SEQRES'] - 1))
c = 0
tot_c = 0
for i in template_ali:
    if c == iter_range:
        break
    if i != '-':
        c += 1
    tot_c += 1
iter_range = tot_c  # account for '-' in the template_ali file to get an accurate starting point
count = abs(min(0, seq_index['structure_begin'] - (seq_index['pseudo_template_SEQRES'] - 1)))
if count != 0:  # because MODELLER expects an alignment to the whole sequence (?)
    append_start = structure_seq[:count]
else:
    append_start = ''
# Counter for the ATOM seq
count2 = abs(min(0, seq_index['structure_begin'] - (seq_index['pseudo_template_SEQRES'] - 1)))
# Counter for the ATOM seq with missing residues set to 'X'
# from where to begin iterating on the pseudo template seq
for i, j in zip(template_ali[iter_range:], query_ali[iter_range:]):
    if count2 == len(structure_seq):  # All ATOM seq residues are done (i.e. no more structure available)
        break
    curr_hetatm = list(structure[0][sys.argv[2]].get_residues())[count]  # get the next residue with structure
    curr_structure = structure_seq[count2]  # get the next residue with structure or a missing 'X'
    if i == '-':  # if the alignment has a gap in the template
        if j != '-':  # do not increment the counters but append the alignment residues
            final_alignment[0].append(i)
            final_alignment[1].append(j)
            continue
        else:  # if the corresponding query also has a gap
            continue
    elif curr_structure == 'X':  # if the structure is missing for this residue
        i = '-'
        count2 += 1  # increment the counter for ATOM seq with missing residues annotated
        if j != '-':
            final_alignment[0].append(i)
            final_alignment[1].append(j)
            continue
        else:
            continue
    elif curr_hetatm.get_resname() in captured_hetatm:  # if a HETATM check if in the list of expected HETATMs
        i = '.'
    else:  # if not a HETATM check if the residue matches in the two ATOM seqs and the template
        temp = abs(min(0, seq_index['structure_begin'] - (seq_index['pseudo_template_SEQRES'] - 1)))
        err_m = 'Sequence mismatch:\n' + structure_seq + '\n' + template_ali
        assert three_to_one(curr_hetatm.get_resname()) == i, err_m
        assert curr_structure == i, 'Sequence mismatch'
    final_alignment[0].append(i)
    final_alignment[1].append(j)
    count += 1
    count2 += 1
if count2 < len(structure_seq):  # MODELLER expects and alignment for the full structure sequence
    append_last = structure_seq[count2:]
else:
    append_last = ''
pdbname = sys.argv[1].split('/')[-1].split('.')[0]
with open(sys.argv[3] + '/alignment_MODELLER_' + sys.argv[4] + '-' + pdbname + '.ali', 'w') as f:
    # output the alignment file (MODELLER ready)
    f.write('>P1;' + pdbname + '\n')
    # appending only a single chain information. For multi-chain collation -> separate script
    f.write('structure' + ':' + pdbname + ':FIRST:' + sys.argv[2] + ':LAST:' + sys.argv[2] + '::::\n')
    f.write(append_start + ''.join(final_alignment[0]) + append_last +  '*\n')
    f.write('>P1;' + sys.argv[4] + '\n')
    f.write('sequence:::::::::\n')
    f.write(('-' * len(append_start)) + ''.join(final_alignment[1]) + ('-' * len(append_last)) + '*')

if os.path.isfile(sys.argv[3] + '/alignment_res_num_table_pseudo.txt'):
    f = open(sys.argv[3] + '/alignment_res_num_table_pseudo.txt', 'a')
else:
    f = open(sys.argv[3] + '/alignment_res_num_table_pseudo.txt', 'w')
    f.write('object\tprot_index_pseudo\tprot_index_MODELLER\tseq_index_pseudo\tseq_index_MODELLER\n')
f.write('query_' + sys.argv[4] + '\t\t\t')
f.write(str(seq_index['pseudo_query_input']) + '\t')
f.write(str(seq_index['pseudo_query_input'] + iter_range) + '\n')
f.write('template_' + sys.argv[1] + '\t')
f.write(str(prot_index['structure_begin'] - seq_index['structure_begin'] +
            seq_index['pseudo_template_SEQRES']) + '\t')
f.write(str(prot_index['structure_begin'] - seq_index['structure_begin'] +
            seq_index['pseudo_template_SEQRES'] + iter_range) + '\t')
f.write(str(seq_index['pseudo_template_SEQRES']) + '\t')
f.write(str(seq_index['pseudo_template_SEQRES'] + iter_range) + '\n')
f.close()
