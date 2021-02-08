import sys
import re

# Satwik Pasani: stochastic13 (29-Jan-2021)

# Adds residue numbering offset and a valid chain name to all of the ATOM entries in the given PDB
# (Necessary for MODELLER outputs before passing to IMP)
# optionally deletes contiguous residue segments >n which do not have a corresponding structure in the template

# Input Arguments (separated by space):
# 1. Name of the PDB
# 2. Chain name (only used if the Chain identifier is not already present.
# Single Chain name is used for all entries otherwise)
# 3. Residue offset (for chain_wise residue offset enter something similar to: "A:12,B:45")
# (For fixing PDB files where a sequence present in the corresponding FASTA is missing from the PDB, this script
# can offset the residue numbers starting from a given value. Run with that option first, and then the actual numbering
# offset. Enter '23;12' to offset all residues >= residue number 23 by 12. This option does not look at the chains.)
# 4. Output file name
# 5, 6: Optional. If present, equates to n and ali_file respectively. Deletes all contiguous segments of length >n
# from the pdb if such a segment is not paired to a corresponding structure in the provided .ali file. It takes the
# sequence corresponding to a "structure" to be the template sequence. Currently ONLY FOR SINGLE CHAIN MODELS

chain_offset = {}
if ':' in sys.argv[3]:  # chain-wise offset
    chain_wise = True  # flag indicating the offset is according to chains
    for i in sys.argv[3].split(','):
        k, v = i.split(':')
        chain_offset[k] = v
else:
    chain_wise = False

cutoff_resnum = -1e7  # arbitrary high value (declaration not needed)
if ';' in sys.argv[3]:  # offset only after a particular residue number
    specific_offset = True
    cutoff_resnum = int(sys.argv[3].split(';')[0])
    sys.argv[3] = sys.argv[3].split(';')[1]  # the actual offset
else:
    specific_offset = False

assert len(sys.argv[2]) == 1, 'Chain Name should be one letter'
with open(sys.argv[1], 'r') as f_in:  # read the PDB
    rd = f_in.read().split('\n')
all_resnums = set()  # later useful for the next section
with open(sys.argv[4], 'w') as f_out:
    for line in rd:
        if len(line) == 0:  # skip empty lines
            continue
        if not (line[:4].strip() in ['ATOM', 'TER']):  # skip all entries apart from ATOM/TER entries
            f_out.write(line + '\n')
            continue
        newline = line[:21]
        chain_name = line[21]  # column 22
        if chain_name == ' ':  # no chain_name specified
            newchain_name = sys.argv[2].upper()
        else:  # chain name already present
            newchain_name = chain_name  # use the old chain name
        newline = newline + newchain_name  # append the chain name
        resnum = int(line[22:26].strip())  # columns 23 - 26
        if specific_offset and (resnum < cutoff_resnum):
            newresnum = str(resnum)  # no offset
        else:  # no specific offset (i.e. all values to be offset), or if beyond the cutoff for specific offset
            if not chain_wise:
                newresnum = str(resnum + int(sys.argv[3]))
            else:
                assert newchain_name in chain_offset, 'Must enter a valid offset for all chains'
                newresnum = str(resnum + int(chain_offset[newchain_name]))
        assert len(newresnum) <= 4, 'Residue number larger than 4 digits'
        newresnum = '{:>4s}'.format(newresnum)  # pad the number appropriately
        newline += newresnum
        newline += line[26:]
        f_out.write(newline + '\n')
        all_resnums.add(newresnum)

# TODO: Implement multi-chain model support
if len(sys.argv) > 5:  # long continuous segment deletion options present
    assert len(sys.argv) == 7, "Invalid number of arguments"
    segment_size = int(sys.argv[5])
    ali_file = sys.argv[6]
    reg = r'\n*([A-Z-.]+\n)*([A-Z-.]+[*]\n*)'  # regex to identify the two sequences
    reg2 = r'^[\s\S]*sequence[\s\S]*structure[\s\S]*$'  # regex to find the order of the sequence and the structure
    reg3 = r'^[\s\S]*structure[\s\S]*sequence[\s\S]*$'  # only one of these two should match
    with open(ali_file, 'r') as ali_file:
        rd = ali_file.read()
    if re.search(reg2, rd) is None:
        assert not (re.search(reg3, rd) is None), 'No Structure/Sequence identifier found in the ali file'
        order = 1  # structure before sequence
    else:
        assert re.search(reg3, rd) is None, '>2 Structure/Sequence identifiers found in the ali file'
        order = 0  # sequence before structure
    # join all the identified regex groups, remove all newlines, strip trailing asterisk
    seq1, seq2 = [''.join(x).strip('*').replace('\n', '') for x in re.findall(reg, rd)]
    if order == 0:
        final_seq = seq2
        compare_seq = seq1  # compare_seq allows identifying leading "-" which need to be removed
    else:
        final_seq = seq1
        compare_seq = seq2
    m = re.match("^-+", compare_seq)
    if m:  # leading "-" need to be removed
        final_seq = final_seq[m.end():]
    reg = '[-]{%d,5000}' % segment_size  # identify (greedily) contiguous "-" segments longer than segment_size
    if re.search(reg, final_seq) is None:  # nothing to be replaced
        print("No segment to be deleted.")
        quit(0)
    to_be_deleted = set()  # all indices to be deleted
    for m in re.finditer(reg, final_seq):  # find all matches of the regex
        start, end = m.span()
        print("Will delete from indices", start, "to", end)
        to_be_deleted |= set([i for i in range(start, end)])
    with open(sys.argv[4], 'r') as f:
        rd = f.read()
        rd = rd.split('\n')
    all_resnums = sorted(list(all_resnums))
    with open(sys.argv[4], 'w') as f:  # overwrite the output file
        for line in rd:
            if len(line) == 0:  # skip empty lines
                continue
            if line[:4] != 'ATOM':  # skip all entries apart from ATOM entries
                f.write(line + '\n')
                continue
            resnum = line[22:26]
            if all_resnums.index(resnum) in to_be_deleted:
                continue  # do not add this to the output
            f.write(line + '\n')
