#!/usr/bin/env python3

import numpy as np

def mask_for_gc4(sequence):
    # check for triplet
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) is not divisible by 3")

    # convert string to a numpy array of ASCII bytes for fast operations
    integer_sequence = np.frombuffer(sequence.upper().encode('ascii'), dtype=np.uint8)
    
    # get the nucleotides for each position
    nuc1 = integer_sequence[0::3]
    nuc2 = integer_sequence[1::3]

    # convert letters to ASCII integers # A 65; T 84; G 71; C 67
    A, C, G, T = ord('A'), ord('C'), ord('G'), ord('T')
    
    # make flag amino acids
    is_val = (nuc1 == G) & (nuc2 == T)
    is_pro = (nuc1 == C) & (nuc2 == C)
    is_thr = (nuc1 == A) & (nuc2 == C)
    is_ala = (nuc1 == G) & (nuc2 == C)
    is_gly = (nuc1 == G) & (nuc2 == G)
    is_leu = (nuc1 == C) & (nuc2 == T)
    is_arg = (nuc1 == C) & (nuc2 == G)
    is_ser = (nuc1 == T) & (nuc2 == C)
    
    # this allows me to evaluate a single codon for any of the following codons
    is_candidate_codon = is_val | is_pro | is_thr | is_ala | is_gly | is_leu | is_arg | is_ser
    
    # make an empty sequence to fill with values of the boolean mask
    full_mask = np.zeros(len(integer_sequence), dtype=bool)
    full_mask[2:len(integer_sequence):3] = is_candidate_codon
    
    return full_mask

print(mask_for_gc4('ATGCCCATT'))


    