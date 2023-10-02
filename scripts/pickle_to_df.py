#!/usr/bin/env python3
import argparse
from collections import Counter

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from tral_score import load_repeatlists

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r", "--repeat_dir", type=str, required=True, help="Path to directory containing pickled RepeatLists"
    )
    parser.add_argument(
        "-o", "--output_file", type=str, required=True, help="Name for output file"
    )

    return parser.parse_args()

def get_consensus_unit(msa: list) -> str:    
    """ Given a list of units representing a multiple sequence
    alignment (from DB), determine the consensus unit at each column of the alignment and 
    return the consensus unit. Insertion columns where the number of 
    insertions >= number of nucleotides are skipped. For each column, most common nt
    is selected as consensus, in case of tie: pick one at random.
    Parameters
    msa (str):  String representation of a multiple sequence alignment, with units
                delimited by commas
    
    Returns  
    consensus_unit (str): 
                A consensus unit derived from the provided msa
    """    
    # Convert msa into list of lists, convert into np.Array() and transpose
    msa_matrix_t = np.array([list(unit) for unit in msa]).transpose()
    consensus_unit = ""

    for msa_col in msa_matrix_t:
        # if half or more of the msa column are gap entries, skip column
        if np.count_nonzero(msa_col == '-') >= 0.5 * len(msa_col):
            continue

        # discard gap entries, get most common nt in column
        msa_col = msa_col[msa_col != '-']
        consensus_unit += Counter(msa_col).most_common(1)[0][0] # most_common() will return e.g. [('A': 5)]

    return consensus_unit

def standardize_unit(unit: str) -> str:
    """ Determine all circular permutations of input unit and its reverse complement, sort 
    alphabetically and return the first option. 
    
    e.g. TTTA -> ['TTTA', 'TAAA', 'TTAT', 'AAAT', 'TATT', 'AATA', 'ATTT', 'ATAA'] -> AAAT
    """
    rev_cmp = str(Seq(unit).reverse_complement())
    permutations = []

    for i in range(len(unit)):
        permutations.append(unit[i:len(unit)] + unit[0:0+i])
        permutations.append(rev_cmp[i:len(unit)] + rev_cmp[0:0+i])

    return sorted(permutations)[0]


def main():
    args = parse_cla()

    df_tandem_repeats = {
        "seq": [],
        "begin": [],
        "l_effective": [],
        "n_effective": [],
        "repeat_region_length": [],
        "score": [],
        "pvalue": [],
        "divergence": [],
        "standardized_consensus_unit": [],
        "msa": [],
    }

    for file_name, repeat_list in load_repeatlists(args.repeat_dir):
        for repeat in repeat_list.repeats:
            consensus_unit = get_consensus_unit(repeat.msa)
            standardized_unit = standardize_unit(consensus_unit)

            df_tandem_repeats["seq"].append(file_name.replace("_refined.pickle", ""))
            df_tandem_repeats["begin"].append(repeat.begin)
            df_tandem_repeats["l_effective"].append(repeat.l_effective)
            df_tandem_repeats["n_effective"].append(repeat.repeat_region_length/ repeat.l_effective )
            df_tandem_repeats["repeat_region_length"].append(repeat.repeat_region_length)
            df_tandem_repeats["score"].append(repeat.d_score["phylo_gap01"])
            df_tandem_repeats["pvalue"].append(repeat.d_pvalue["phylo_gap01"])
            df_tandem_repeats["divergence"].append(repeat.d_divergence["phylo_gap01"])
            df_tandem_repeats["standardized_consensus_unit"].append(standardized_unit)
            df_tandem_repeats["msa"].append(", ".join(repeat.msa))
    
    df_tandem_repeats = pd.DataFrame(df_tandem_repeats)
    df_tandem_repeats.to_csv(args.output_file, index=False, sep="\t")

if __name__ == "__main__":
    main()
