import vdj.io
import numpy as np




def nucleotide_idx():
    """
    Returns the dictionary linking base identity to an integer
    """
    return {'A':0, 'C':1, 'G':2, 'T':3}

def endogenous_seqs():
    """
    Returns a dictionary of the sequence identity for the reference 12RSS
    sequence and useful mutations. 
    """
    conv = nucleotide_idx()
    _seqs = {'reference': 'CACAGTGCTACAGACTGGAACAAAAACC',
            'V1-135':  'CACAGTGATTCAGACCCGAACAAAAACT',
            'V9-120': 'CACAGTGATACAAATCATAACATAAACC',
            'V10-96': 'CACAATGATATAAGTCATAACATAAACC',
            'V10-95': 'CACAATGATATAAGTCATAACATAAACC',
            'V19-93': 'CACAGTGATACAAATCATAACAAAAACC',
            'V4-55': 'CACAGTGATACAGACTGGAACAAAAACC',
            'V5-43': 'CACAGTGATGCAGACCATAGCAAAAATC',
            'V6-15': 'CACAGTACTTCAGCCTCCTACATAAACC',
            'DFL-161': 'CACAGTGCTATATCCATCAGCAAAAACC',
            'DFL-1613': 'CACAGTAGTAGATCCCTTCACAAAAAGC'}

    seqs = {m: [seq, np.array([conv[a] for a in seq])] for m, seq in _seqs.items()}
    return seqs

# Load the raw sequences and nt conversion id's
conversion = nucleotide_idx()
seqs = endogenous_seqs()
ref = seqs['reference'][0]


mut_id = '12SpacG10CG11A'

# Determine the region which is mutated
if 'hept' in mut_id.lower():
    seq = ref[:7]
    region = 'hept'
elif 'spac' in mut_id.lower():
    seq = ref[7:19]
    region =  'spac'
elif 'non' in mut_id.lower():
    seq = ref[19:]
    region = 'non'
else:
    try:
        new_seq = seqs[mut_id]
        new_seq_idx = np.array([conversion[a] for a in new_seq])
    except:
        raise ValueError("Mutation ID not in proper format of 12(region)(old)(pos)(new) or name not recognized.")


# Force the string to uppercae
loc = mut_id.lower().split(region)[-1].upper()

# Get the indices of base positions
base_id = np.array([i for i, b in enumerate(list(loc)) if b in 'ATCG']).astype(int)

# Split the indices into pairs
if (len(base_id) / 2) == 1:
    muts = [loc]
else:
    inds = [(start, end) for start, end in zip(base_id[::2], base_id[1::2])]
    muts = [loc[ind[0]:ind[1] + 1] for ind in inds]

new_region = list(seq)
for m in muts:
    pos = int(m[1:-1])
    if seq[pos - 1] != m[0].upper():
        raise ValueError(f'Base at position {m[1]} is not {m[0].upper()}! Double check the sequence.')
    else:
        new_region[int(m[1]) - 1] = m[-1]


# Replace the region and compute the integer representation
new_seq = ref.replace(seq, ''.join(new_region).upper())
new_seq_idx = np.array([conversion[a] for a in new_seq]) 
