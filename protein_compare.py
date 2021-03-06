"""This module can read a .txt file of protein sequences and return how similar
they are to each other as well as search for potential structural motifs, such
as Schellman C-caps. This comparison has been designed to analyse the variable
binding loop regions in antibody-like proteins, hence the usage of the term
"loops" in the code. Proteins are allowed to have different numbers of loops
in the input file, but only those loops which exist in all proteins can be
compared.
"""

from re import search, findall, IGNORECASE
from sys import version_info
from decimal import getcontext, Decimal
from itertools import product
from collections import OrderedDict
from collections.abc import Iterable
import matplotlib.pyplot as plt
import pandas as pd

def full_prot_comparison(input_data, **kwargs):
    """A top-level function that users should be interacting with. It compares
    each protein of interest against every protein in the input file.
    Frameshift can be introduced into the comparisons.

    Args:
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a pandas DataFrame of all the protein
            sequences for comparison

    Kwargs:
        shift (tuple): how much frame shift do you want in the comparison.
        proteins_of_interest (list): the protein(s) against which you want to
            run the comparisons. By default all proteins are compared.
        loops (list): which loops do you want compare. Default is just loop 1.

    Returns:
        tuple: a three-member tuple, comprising of a: string, list, and
            how_to_compare() generator
    """

    def how_to_compare():
        """A generator function that specifies the details of how the proteins
        ought to be compared against each other and how the results should be
        labelled.

        Yields:
            tuple: a two-member tuple, the first member being a string label
                for the y-axis, and the other a generator yielding similarity
                scores of comparisons
        """
        for prot_of_interest in proteins_of_interest:
            comparison_scores = compare_proteins(loops, prot_seq, prot_of_interest, shift)
            yield prot_of_interest, comparison_scores

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.index
    proteins_of_interest = kwargs.get("proteins_of_interest", all_prot_list)
    shift = kwargs.get("shift", ("", "", ""))
    loops = kwargs.get("loops", ["loop_1"])
    title = "{} {} shift{}".format(data_source, loops, shift)
    comparison_results = how_to_compare()
    return title, all_prot_list, comparison_results


def one_prot_all_shifts(prot_of_interest, input_data):
    """A top-level function that users should be interacting with as an
    alternative to full_prot_comparison(). This takes a single protein of
    interest, subjects it to all possible shifts in all loops (individually),
    and compares it against all other proteins.

    Args:
        prot_of_interest (str): the name of the protein of interest
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a pandas DataFrame of all the protein
            sequences for comparison

    Returns:
        tuple: a three-member tuple, comprising of a: string, list, and
            how_to_compare() generator
    """

    def how_to_compare():
        """A generator function that specifies the details of how the proteins
        ought to be compared against each other and how the results should be
        labelled.

        Yields:
            tuple: a two-member tuple, the first member being a string label
                for the y-axis, and the other a generator yielding similarity
                scores of comparisons
        """
        loop_number = len(prot_seq.loc[prot_of_interest].dropna())
        loops = ["loop_" + str(i) for i in range(1, loop_number + 1)]
        all_shifts = product(loops, ("first", "second"), range(8))
        for shift in all_shifts:
            comparison_scores = compare_proteins(loops, prot_seq, prot_of_interest, shift)
            yield str(shift), comparison_scores

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.index
    title = "{} protein {} all shifts".format(data_source, prot_of_interest)
    comparison_results = how_to_compare()
    return title, all_prot_list, comparison_results


def all_prots_all_shifts(input_data):
    """A top-level function that users should be interacting with. This takes
    all proteins in a file, subjects them to all possible shifts in all loops
    (individually), and compares them against all other proteins. Only the
    highest comparison score from each protein vs protein matchup is shown,
    compressing the data from 3D to 2D.

    Args:
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a pandas DataFrame of all the protein
            sequences for comparison

    Returns:
        tuple: a three-member tuple, comprising of a: string, list, and
            how_to_compare() generator
    """

    def how_to_compare():
        """A generator function that specifies the details of how the proteins
        ought to be compared against each other and how the results should be
        labelled.

        Yields:
            tuple: a two-member tuple, the first member being a string label
                for the y-axis, and the other a generator yielding similarity
                scores of comparisons
        """
        for prot in prot_seq.index:
            comparison_function = one_prot_all_shifts(prot, (input_data, prot_seq))
            # The title and y axis labels that come out of one_prot_all_shifts are ignored.
            _, result_df = turn_into_dataframe(comparison_function)
            highest_scores_df = result_df.max(axis=0)
            # Highest_scores is a generator object, for the sake of consistency with all the other
            # top-level functions.
            highest_scores = (score for score in highest_scores_df)
            yield prot, highest_scores

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.index
    title = "highest scores from all_prots_all_shifts, {}".format(data_source)
    highest_comparison_results = how_to_compare()
    return title, all_prot_list, highest_comparison_results


def amino_acid_properties_matrix():
    """A top-level function that users should be interacting with. It creates a
    pairwise amino acid comparison matrix that can be directly compared with
    the Blosum62 substitution matrix.

    Returns:
        tuple: a three-member tuple, comprising of a: string, tuple, and
            how_to_compare() generator
    """

    def how_to_compare():
        """A generator function that specifies the details of how the amino
        acids ought to be compared against each other and how the results
        should be labelled.

        Yields:
            tuple: a two-member tuple, the first member being a string label
                for the y-axis, and the other a generator yielding similarity
                scores of comparisons
        """
        comparison_func = lambda aa_2: round_floats(compare_amino_acids(aa_1, aa_2))
        for aa_1 in all_amino_acids:
            comparison_scores = (comparison_func(aa_2) for aa_2 in all_amino_acids)
            yield aa_1, comparison_scores

    all_amino_acids = tuple("ACDEFGHIKLMNQPRSTVWY")  # A tuple of the individual amino acids.
    title = "amino acid comparison matrix"
    comparison_results = how_to_compare()
    return title, all_amino_acids, comparison_results


def round_floats(number):
    """A function which converts float values of comparison scores into floats
    with no more than two decimal figures. No precision is lost this way - the
    point is to convert numbers like 1.7499999999 into 1.75.

    Arguments:
        number (float): the value of a comparison score

    Returns:
        float: value of a comparison score with two decimal figures at most
    """
    getcontext()
    two_decimal_places = Decimal('0.01')
    rounded_number = Decimal(number).quantize(two_decimal_places)
    return float(rounded_number)


def read_data(input_file):
    """This is used to convert a chosen file .txt with sequence data into a
    pandas DataFrame that the rest of this script can work with. The file needs
    to have data separated by tabs or spaces, with each line corresponding to a
    single protein. The peptide comparison algorithm is case-insensitive, but
    requires usage of the single-letter amino acid code.

    Args:
        input_file (str): path to the input file with the protein sequences

    Returns:
        DataFrame: pandas DataFrame object containing all the data from the
            chosen file
    """
    # In python versions 3.7+ the normal dictionary is already ordered, but in lower versions it is
    # necessary to use the OrderedDict object.
    prot_dict = dict() if version_info >= (3, 7) else OrderedDict()

    with open(input_file, "r") as data:
        for line in data:
            words = line.split()
            # A safeguard that ignores lines composed entirely out of whitespace characters.
            if words:
                protein_name, peptides = words[0], words[1:]
                prot_dict[protein_name] = peptides

    prot_seq = pd.DataFrame.from_dict(prot_dict, dtype=str, orient="index")
    loop_number = prot_seq.shape[1]
    prot_seq.columns = ["loop_" + str(i) for i in range(1, loop_number + 1)]
    return prot_seq


def compare_proteins(loops, prot_seq, prot_of_interest, shift):
    """A generator function comparing of one protein of interest against all
    other proteins.

    Args:
        loops (list): which loops do you want compare
        prot_seq (DataFrame): pandas DataFrame object containing all the
            protein peptide sequences
        prot_of_interest (str): the protein against which others are analysed
        shift (tuple): how much frame shift do you want in the comparison. The
            first item in the tuple specifies which loop(s) the shift should
            occur in.

    Yields:
        float: a numerical score of each protein's similarity to the protein of
            interest
    """
    for prot in prot_seq.index:
        comparison_score = 0
        for loop in loops:
            peptide_1, peptide_2 = prot_seq.at[prot_of_interest, loop], prot_seq.at[prot, loop]
            # If loop_2 is short, the score from the loop_1 comparison becomes twice as important,
            # so it is copied.
            if len(peptide_1) < 4 or len(peptide_2) < 4:
                comparison_score += comparison_score
            else:
                shift_peptide, shift_index = (shift[1], shift[2]) if loop in shift[0] else ("", "")
                comparison_score += compare_peptides(peptide_1, peptide_2, shift_peptide,
                                                     shift_index)
        yield round_floats(comparison_score)


def compare_peptides(peptide_1, peptide_2, shift_peptide, shift_index):
    """This compares the sequences of two peptides for similarities. If
    compare_peptides() was called by compare_proteins() then peptide_1 comes
    from the protein of interest and peptide_2 from the comparison protein.

    Args:
        peptide_1 (str): first of the peptides being compared, written in
            single-letter code
        peptide_2 (str): second of the peptides being compared, written in
            single-letter code
        shift_peptide (str): which peptide to place the X amino acid in
        shift_index (int): which index to place the X amino acid in

    Returns:
        float: a numerical score of how similar the peptides are to each other
    """
    # This code block is verbose, but optimised for speed.
    if shift_peptide == "first":
        peptide_1_iter = iterate_peptide(peptide_1, shift_index)
    else:
        peptide_1_iter = iter(peptide_1)
    if shift_peptide == "second":
        peptide_2_iter = iterate_peptide(peptide_2, shift_index)
    else:
        peptide_2_iter = iter(peptide_2)

    # If the peptides are of different lengths, this zip() function will truncate the longer peptide
    # down to the length of the shorter one.
    amino_acid_pairs = zip(peptide_1_iter, peptide_2_iter)
    similarity_scores = (compare_amino_acids(*aa_pair) for aa_pair in amino_acid_pairs)
    return sum(similarity_scores)


def iterate_peptide(peptide, shift_index):
    """A generator function, which yields a peptide one amino-acid at a time.
    It inserts an additional X amino acid into the peptide at a given position,
    which results in frameshift.

    Args:
        peptide (str): the peptide being split into amino acids
        shift_index (int): which index in the peptide an "X" amino acid should
            be inserted in

    Yields:
        str: one amino-acid written in the single-letter code
    """
    for index, character in enumerate(peptide):
        if index == shift_index:
            yield "X"
        yield character


def compare_amino_acids(aa_1, aa_2):
    """This compares two amino acids for similarities in terms of size,
    hydrophobicity, and miscallaneous qualities.

    Args:
        aa_1 (str): the single-letter code of the first amino-acid
        aa_2 (str): the single-letter code of the second amino-acid

    Returns:
        float: numerical score of how similar the amino acids are to each other
    """
    # Making sure that the strings are all upper case.
    aa_1, aa_2 = aa_1.upper(), aa_2.upper()
    if aa_1 == aa_2 != "X":
        # G, P, and C amino acids are more important to preserve, hence they get the higher score.
        misc = 3 if aa_1 in "GPC" else 2
    else:
        # This is the value both for non-matching amino acids and matching "X" amino acids.
        misc = 0
        for category in MISC_SCORES:
            if aa_1 in category[0] and aa_2 in category[0]:
                misc = category[1]
                break  # To not waste time.

    hydro = MAX_HYDRO - abs(HYDRO_SCORES[aa_1] - HYDRO_SCORES[aa_2])
    size = MAX_SIZE - abs(SIZE_SCORES[aa_1] - SIZE_SCORES[aa_2])
    return misc * MISC_WEIGHT + hydro * HYDRO_WEIGHT + size * SIZE_WEIGHT


def fragmentate_dictionary(old_dictionary):
    """This splits dictionary entries like "FW": 0 into "F":0 and "W": 0.

    Args:
        old_dictionary (dict): a dictionary containing either hydrophobicity or
            size scores for amino acids

    Returns:
        dict: a less-readable but more computationally efficient dictionary,
            compared to the original one
    """
    new_dictionary = {}
    for key, value in old_dictionary.items():
        for character in key:
            new_dictionary[character] = value
    return new_dictionary


HYDRO_SCORES = fragmentate_dictionary({"FW": 0, "ILMY": 0.75, "V": 1.5, "AP": 2.25, "XG": 3,
                                       "STC": 4, "NQH": 5, "DERK": 6})
MAX_HYDRO = 6
HYDRO_WEIGHT = 1.0

SIZE_SCORES = fragmentate_dictionary({"GASC": 0, "VTPDN": 1, "X": 1.5, "EQILMHFK": 2, "WYR": 3})
MAX_SIZE = 3
SIZE_WEIGHT = 1.2

# These are respectively: positively charged, negatively charged, Asx/Glx, beta-branched, and
# neutral aromatic. "YF" needs to be before "HYWF".
MISC_SCORES = (("RHK", 1), ("DE", 1), ("ND", 1.5), ("QE", 1.5), ("TVLI", 1), ("YF", 1.5),
               ("HYWF", 1))
MISC_WEIGHT = 1.8


def write_csv_file(output_file, comparison_function):
    """A top-level function that turns the outputs of other top-level functions
    into .csv (comma-separated values) files that can be opened with Excel.

    Args:
        comparison_function: which function's results should get converted into
            an csv file
        output_file (str): path to the file that the data gets put out into
    """
    title, result_df = turn_into_dataframe(comparison_function)
    with open(output_file + ".csv", "a+") as csv_file:
        csv_file.write(title)
    result_df.to_csv(output_file + ".csv", mode="a+")


def quick_analysis(input_file, output_file):
    """A function to make analysis quicker. Runs a protein comparison in three
    frames (0, +1 for loop 1, and +1 for loop 2) and writes the results into a
    single .txt file.

    Args:
        input_file (str): path to the input file with the protein sequences
        output_file (str): path to the file that the data gets put out into
    """
    input_data = (input_file, read_data(input_file))
    loops = ["loop_1", "loop_2"]
    analysis_func = lambda shift: write_csv_file(output_file, full_prot_comparison(input_data,
            loops=loops, shift=shift))
    for shift in (("", "", ""),
                  (["loop_1"], "first", 1),
                  (["loop_1"], "second", 1),
                  (["loop_2"], "first", 1),
                  (["loop_2"], "second", 1)):
        analysis_func(shift)


def unpack_generators(data, unwanted_types=None):
    """A utility function which recursively searches for generator objects
    within an iterable and turns them into tuples. Optionally, it also confirms
    that the iterable contains no un-expected data object types.

    Args:
        data (any): often an iterable containing generator objects, but can be
            anything
        unwanted_types (list): a list of all the object types that should not
            occur within data. None by default, meaning no screening for
            unwanted object types is done.

    Returns:
        tuple or any: a tuple containing only tuples, strings, and
            non-iterables. If the data argument was a string or not an
            iterable, it is returned un-changed.
    """
    if unwanted_types:
        for unwanted_type in unwanted_types:
            if isinstance(data, unwanted_type):
                raise TypeError("unwanted object type {} found: {}".format(unwanted_type, data))

    if not isinstance(data, str) and isinstance(data, Iterable):
        # As a side-effect this converts all lists, dictionaries, sets, etc. into tuples as well.
        return tuple(unpack_generators(item, unwanted_types) for item in data)
    else:
        return data


def turn_into_dataframe(comparison_function):
    """A top-level function that creates a pandas DataFrame object out of the
    output of other top-level functions.

    Args:
        comparison_function: which function's results should get converted into
            a DataFrame

    Returns:
        tuple: a two-member tuple, the first member being a descriptive string
            and the other being a DataFrame object containing all the
            comparison scores
    """
    title, x_axis_label, comparison_results = comparison_function
    # Unpacking the results is necessary only for the amino_acid_properties_matrix function - all
    # the others are fine without this line.
    unpacked_result = unpack_generators(comparison_results)
    # In python versions 3.7+ the normal dictionary is already ordered, but in lower versions it is
    # necessary to use the OrderedDict object.
    result_dict = dict(unpacked_result) if version_info >= (3, 7) else OrderedDict(unpacked_result)
    result_df = pd.DataFrame.from_dict(result_dict, columns=x_axis_label, orient="index")
    return title, result_df


def find_structural_motifs(peptide):
    """This searches for structural motifs, like alpha helix N-caps and salt
    bridges, in the provided peptide.

    Args:
        peptide (str): the peptide being analysed

    Returns:
        list: a list of all the motifs found within the peptide
    """

    def scan(motif_dictionary):
        """A generator function that scans a peptide for motifs found in the
        specified motif dictionary.

        Args:
            motif_dictionary (dict): a dictionary of regex search patterns as
                keys and the corresponding structural motifs as values.
                Sometimes the values contain further functions to specify the
                particular sub-type of motif that was found.

        Yields:
            str: names of motifs found within the peptide
        """
        for motif_pattern, motif_name in motif_dictionary.items():
            if search(motif_pattern, peptide, IGNORECASE):
                yield motif_name

    stringency = lambda: "Narrow" if search("[LAM]..[HKL]G[IVK]", peptide, IGNORECASE) else "Broad"
    box_ncap_variants = {"[MLIFV][STN]..[EQ][MFIFV]": " with hydrophobic interactions",
                         "[STN]P.[EQ]": " with a proline"}
    big_box_ncap_variants = {"[MLIFV][STN]...[EQ][MFIFV]": " with hydrophobic interactions",
                             "[STN]P..[EQ]": " with a proline"}

    structural_motifs = {"[STN]..[EQ]": "Box N-cap" + "".join(scan(box_ncap_variants)),
                         "[STN]...[EQ]": "Big Box N-cap" + "".join(scan(big_box_ncap_variants)),
                         "[RK]...[E]": "i-->i+4 Salt bridge",
                         "[E]...[RK]": "i-->i+4 Salt bridge",
                         "[RK]..[E]": "i-->i+3 Salt bridge",
                         "[E]..[R]": "i-->i+3 Salt bridge",
                         "[IVKMLYFWA]...G[GSTNEDQ][IVKMLYFWA]": "AlphaL C-cap",
                         "[IVKMLYFWA]...G[IVKMLYFWA]": stringency() + " Shellman C-cap"}

    found_motifs = list(scan(structural_motifs))
    if len(findall("[TVI]", peptide, IGNORECASE)) > 1:
        found_motifs.append("Multiple beta-branched residues")
    return found_motifs


def prot_structural_motifs(input_file):
    """A top-level function that users should be interacting with. It displays
    the potential structural motifs within a set of proteins.

    Args:
        input_file (str): path to the input file with the protein sequences
    """
    prot_seq = read_data(input_file)
    for prot in prot_seq.index:
        print(prot, prot_seq.at[prot, "loop_1"], " ", prot_seq.at[prot, "loop_2"])
        print(find_structural_motifs(prot_seq.at[prot, "loop_1"]),
              find_structural_motifs(prot_seq.at[prot, "loop_2"]))

#pylint:disable=invalid-name
def heatmap(comparison_function, **kwargs):
    """A top-level function for creating heatmaps from protein comparisons.

    Args:
        comparison_function: which function's results should get converted into
            a heatmap plot. Can accept both three-member tuples (outputs from
            functions like one_prot_all_shifts) and two-member tuples (output
            of turn_into_dataframe function).

    Kwargs:
        self_comparisons (bool): whether to display the high similarity scores
            of proteins/amino-acids being compared against themselves. True by
            default.
        map_colours (str): what colour scheme should be used for the heatmap.
            Coolwarm by default.
    """
    self_comparisons = kwargs.get("self_comparisons", True)
    map_colours = kwargs.get("map_colours", "coolwarm")
    if len(comparison_function) == 2:
        title, result_df = comparison_function
    else:
        title, result_df = turn_into_dataframe(comparison_function)

    if not self_comparisons:
        try:
            # Returns anything of any length between "protein " and " all shifts".
            column_name = search(r"protein (.+) all shifts$", title).group(1)
            for index_name in result_df.index:
                result_df.at[index_name, column_name] = None

        except AttributeError:
            # Must use result_df.index instead of result_df.columns, as sometimes the dataframes's
            # shape is rectangular.
            for index_name in result_df.index:
                result_df.at[index_name, index_name] = None

    # Creating the heatmap and colorbar.
    plt.imshow(result_df, cmap=map_colours)
    plt.colorbar()
    # Adding axis units, axis labels, and heatmap title.
    plt.xticks(range(len(result_df.columns)), result_df.columns)
    plt.xlabel("all proteins")
    plt.yticks(range(len(result_df.index)), result_df.index)
    plt.ylabel("variable")
    plt.title(title)
    # The below line is necessary due to a bug in matplotlib - will be removed in the future.
    plt.ylim(len(result_df.columns)-0.5, -0.5)
    plt.show()

###############################################################################

# Here are some examples of how to use the functions found in this module.
if __name__ == "__main__":
    name, data = "example_data.txt", read_data("example_data.txt")

    prot_structural_motifs(name)

    write_csv_file("amino_acid_pairwise_matrix", amino_acid_properties_matrix())

    heatmap(one_prot_all_shifts("prot_two", [name, data]), self_comparisons=False)

    heatmap(full_prot_comparison([name, data], loops=["loop_1", "loop_2"],
            shift=(["loop_1"], "second", 3)), self_comparisons=False)

    title, result_df = turn_into_dataframe(all_prots_all_shifts([name, data]))
    print("\n", result_df)
    heatmap([title, result_df])
