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
from itertools import count, product
from collections import OrderedDict
from collections.abc import Iterable
import matplotlib.pyplot as plt
import numpy as np

def full_prot_comparison(input_data, **kwargs):
    """A top-level function that users should be interacting with. It compares
    each protein of interest against every protein in the input file.
    Frameshift can be introduced into the comparisons.

    Args:
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a dictionary of all the protein peptide
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
            tuple: a two-member tuple, the first member being a generator
                yielding similarity scores of comparisons, and the other being
                the string label for the y-axis
        """
        for prot_of_interest in proteins_of_interest:
            comparison_scores = compare_proteins(loops, prot_seq, prot_of_interest, shift)
            yield (comparison_scores, prot_of_interest)

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.keys()
    proteins_of_interest = kwargs.get("proteins_of_interest", all_prot_list)
    shift = kwargs.get("shift", ("", "", ""))
    loops = kwargs.get("loops", ["loop_1"])
    header = "{} {} shift{}".format(data_source, loops, shift)
    comparison_results = how_to_compare()
    return header, all_prot_list, comparison_results


def one_prot_all_shifts(prot_of_interest, input_data):
    """A top-level function that users should be interacting with as an
    alternative to full_prot_comparison(). This takes a single protein of
    interest, subjects it to all possible shifts in all loops (individually),
    and compares it against all other proteins.

    Args:
        prot_of_interest (str): the name of the protein of interest
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a dictionary of all the protein peptide
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
            tuple: a two-member tuple, the first member being a generator
                yielding similarity scores of comparisons, and the other being
                the string label for the y-axis
        """
        loop_number = len(prot_seq[prot_of_interest])
        loops = ["loop_" + str(i) for i in range(1, loop_number + 1)]
        all_shifts = product(loops, ("first", "second"), range(8))
        for shift in all_shifts:
            comparison_scores = compare_proteins(loops, prot_seq, prot_of_interest, shift)
            yield comparison_scores, str(shift)

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.keys()
    header = "{} protein {} all shifts".format(data_source, prot_of_interest)
    comparison_results = how_to_compare()
    return header, all_prot_list, comparison_results


def all_prots_all_shifts(input_data):
    """A top-level function that users should be interacting with. This takes
    all proteins in a file, subjects them to all possible shifts in all loops
    (individually), and compares them against all other proteins. Only the
    highest comparison score from each protein vs protein matchup is shown,
    compressing the data from 3D to 2D.

    Args:
        input_data (tuple): a two-member tuple containing: a string describing
            the source of the data, and a dictionary of all the protein peptide
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
            tuple: a two-member tuple, the first member being a generator
                yielding similarity scores of comparisons, and the other being
                the string label for the y-axis
        """
        for prot in prot_seq:
            comparison_results = one_prot_all_shifts(prot, (input_data, prot_seq))
            # The header and y axis labels that come out of one_prot_all_shifts are ignored.
            highest_scores = select_highest_scores(comparison_results)
            yield highest_scores, prot

    data_source, prot_seq = input_data
    all_prot_list = prot_seq.keys()
    header = "highest scores from all_prots_all_shifts, {}".format(data_source)
    highest_comparison_results = how_to_compare()
    return header, all_prot_list, highest_comparison_results


def select_highest_scores(comparison_results):
    """A generator function that compresses the output of a top-level function,
    like one_prot_all_shifts. It yields the highest similarity score from each
    column.

    Arguments:
        comparison_results (tuple): a three-member tuple, comprising of a
            string, list, and a how_to_compare() generator

    Yields:
        str: numerical value of the highest similarity score from one column
    """
    unpacked_data = unpack_generators(comparison_results[2])
    for i in count():
        try:
            highest_score = 0
            for scores, _ in unpacked_data:
                if int(scores[i]) > highest_score:
                    highest_score = int(scores[i])
            yield str(highest_score)
        except IndexError:
            break


def amino_acid_properties_matrix():
    """A top-level function that users should be interacting with. It creates a
    pairwise amino acid substitution matrix that can be directly compared with
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
            tuple: a two-member tuple, the first member being a generator
                yielding similarity scores of comparisons, and the other being
                the string label for the y-axis
        """
        comparison_func = lambda aa_2: str(round(compare_amino_acids(aa_1, aa_2)))
        for aa_1 in all_amino_acids:
            comparison_scores = (comparison_func(aa_2) for aa_2 in all_amino_acids)
            yield comparison_scores, aa_1

    all_amino_acids = tuple("ACDEFGHIKLMNQPRSTVWY")  # A tuple of the individual amino acids.
    header = "amino acid comparison matrix"
    comparison_results = how_to_compare()
    return header, all_amino_acids, comparison_results


def read_data(input_file):
    """This is used to convert a chosen file .txt with sequence data into a
    dictionary that the rest of this script can work with. The file needs to
    have data separated by tabs or spaces, with each line corresponding to a
    single protein. The peptide comparison algorithm is case-insensitive, but
    requires usage of the single-letter amino acid code.

    Args:
        input_file (str): path to the input file with the protein sequences

    Returns:
        dict: a dictionary containing all the data from the chosen file
    """
    # In python versions 3.7+ the normal dictionary is already ordered, but in lower versions
    # it is necessary to use the OrderedDict object.
    prot_seq = dict() if version_info >= (3, 7) else OrderedDict()

    with open(input_file, "r") as data:
        for line in data:
            words = line.split()
            protein_name, peptides = words[0], words[1:]
            # This generator expression must be inside the for loop.
            loops = ("loop_" + str(i) for i in count(1))
            # Create a sub-dictionary of peptides as values with keys like "loop_1".
            prot_seq[protein_name] = dict(zip(loops, peptides))
    return prot_seq


def compare_proteins(loops, prot_seq, prot_of_interest, shift):
    """A generator function comparing of one protein of interest against all
    other proteins.

    Args:
        loops (list): which loops do you want compare
        prot_seq (dict): a dictionary of all the protein peptide sequences
        prot_of_interest (str): the protein against which others are analysed
        shift (tuple): how much frame shift do you want in the comparison. The
            first item in the tuple specifies which loop(s) the shift should
            occur in.

    Yields:
        str: a numerical score of each protein's similarity to the protein of
            interest
    """
    for prot in prot_seq:
        comparison_score = 0
        for loop in loops:
            peptide_1, peptide_2 = prot_seq[prot_of_interest][loop], prot_seq[prot][loop]
            # If loop_2 is short, the score from the loop_1 comparison becomes twice as important,
            # so it is copied.
            if len(peptide_1) < 4 or len(peptide_2) < 4:
                comparison_score += comparison_score
            else:
                shift_peptide, shift_index = (shift[1], shift[2]) if loop in shift[0] else ("", "")
                comparison_score += compare_peptides(peptide_1, peptide_2, shift_peptide,
                                                     shift_index)
        yield str(round(comparison_score))


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
        int: a numerical score of how similar the peptides are to each other
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
        int: a numerical score of how similar the amino acids are to each other
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


def write_excel_file(output_file, comparison_function):
    """A top-level function that turns the outputs of other top-level functions
    into .txt files, formatted in such as way as to then easily turn them into
    excel files.

    Args:
        comparison_function: which function's results should get converted into
            an excel file
        output_file (str): path to the file that the data gets put out into
    """

    def data_stream_generator():
        """A generator function which orders the analysis data into strings,
        each string corresponding to a single line in the new file.
        """
        yield header + "\n"
        yield "\t{}\n".format("\t".join(x_axis_label))
        for comparison_scores, y_axis_label in comparison_results:
            yield "{}\t{}\n".format(y_axis_label, "\t".join(comparison_scores))
        yield "\n"

    header, x_axis_label, comparison_results = comparison_function
    with open(output_file + ".txt", "a+") as new_file:
        new_file.writelines(data_stream_generator())


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
    analysis_func = lambda shift: write_excel_file(output_file, full_prot_comparison(input_data,
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
    for prot in prot_seq:
        print(prot, prot_seq[prot]["loop_1"], " ", prot_seq[prot]["loop_2"])
        print(find_structural_motifs(prot_seq[prot]["loop_1"]),
              find_structural_motifs(prot_seq[prot]["loop_2"]))

#pylint:disable=invalid-name
def heatmap(comparison_function, **kwargs):
    """A top-level function for creating heatmaps from protein comparisons.

    Args:
        comparison_function: which function's results should get converted into
            a heatmap plot

    Kwargs:
        self_comparisons (str): used to remove the high similarity scores of
            proteins/amino-acids being compared against themselves. For the
            amino_acid_matrix and full_prot_comparison use "diagonal", and for
            one_prot_all_shifts use "vertical".
        map_colours (str): what colour scheme should be used for the heatmap.
            Coolwarm by default.
    """
    self_comparisons = kwargs.get("self_comparisons", None)
    map_colours = kwargs.get("map_colours", "coolwarm")
    all_data = unpack_generators(comparison_function)
    title, x_labels = all_data[0], all_data[1]
    # Splitting one list of two-member tuples into two tuples by using the *
    scores, y_labels = zip(*all_data[2])
    # A 2D numpy array, needs to be an array of floats for image formation.
    score_array = np.array(scores, float)


    y_indexes = range(len(y_labels))
    if self_comparisons == "diagonal":
        for index in y_indexes:
            # Setting cells to None raises a warning at runtime, but the function still works fine.
            score_array[index][index] = None

    elif self_comparisons == "vertical":
        # Returns anything of any length between "protein " and " all shifts".
        prot_of_interest = search(r"protein (.+) all shifts$", title).group(1)
        x_index = x_labels.index(prot_of_interest)
        for y_index in y_indexes:
            score_array[y_index][x_index] = None

    elif self_comparisons is not None:
        raise ValueError("Invalid value for keyword argument self_comparisons. Choose either" +
            "'vertical', 'diagonal', or None.")

    # Creating the heatmap and colorbar.
    fig, ax = plt.figure(), plt.subplot()
    im = ax.imshow(score_array, cmap=map_colours)
    fig.colorbar(im, ax=ax)
    # Adding axis labels, axis units, and heatmap title.
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_xticklabels(x_labels)
    ax.set_yticklabels(y_labels)
    plt.xlabel("protein of comparison")
    plt.ylabel("variable")
    ax.set_title(title)
    # The below line is necessary due to a bug in matplotlib - will be removed in the future.
    ax.set_ylim(len(y_labels)-0.5, -0.5)

    plt.show()

###############################################################################

# Here are some examples of how to use the functions found in this module.
if __name__ == "__main__":
    exmpl = "example_data.txt"

    prot_structural_motifs(exmpl)

    write_excel_file("amino_acid_pairwise_matrix", amino_acid_properties_matrix())

    heatmap(one_prot_all_shifts("prot_two", (exmpl, read_data(exmpl))), self_comparisons="vertical")

    heatmap(full_prot_comparison((exmpl, read_data(exmpl)), loops=["loop_1", "loop_2"],
            shift=(["loop_1"], "second", 3)), self_comparisons="diagonal")

    heatmap(all_prots_all_shifts((exmpl, read_data(exmpl))), self_comparisons="diagonal")
