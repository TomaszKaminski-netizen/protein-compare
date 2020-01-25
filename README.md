# ProtCompare

ProtCompare is a simple module for comparing proteins by their primary amino
acid sequence, as to determine how similar they are to each other.

## Setup

Simply download the protein_compare.py and example_data.txt into the same folder.
As the name suggests, example_data.txt contains some protein sequences formatted
in a manner which protein_compare.py functions can process.

## Usage

There are four top-level functions that users should be interacting with:
    one_prot_all_shifts
    full_peptide_comparison
    amino_acid_properties_matrix
    prot_structural_motifs

The results of the first three of these functions can be passed on onto higher-
order functions which either populate a .txt file with the results or turn them
into a matplotlib heatmap.

## License
[MIT](https://choosealicense.com/licenses/mit/)
