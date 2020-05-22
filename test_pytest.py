"""This module is used for unit testing the protein_compare module."""

from pytest import raises
import pandas as pd
from protein_compare import one_prot_all_shifts, unpack_generators, full_prot_comparison, \
    compare_peptides, amino_acid_properties_matrix, find_structural_motifs, all_prots_all_shifts
# pylint:disable=missing-function-docstring

def test_unpack_generators():
    # Checking that unpack_generators correctly raises exceptions.
    bad_test_data = ("aaa", (i for i in range(6)), {"b": 0, "c": 1, "d": 2})
    with raises(TypeError):
        unpack_generators(bad_test_data, [dict])

    valid_test_data = ("aaa", (i for i in range(6)), ("bbb", (i for i in range(4))))
    expected_result = ("aaa", (0, 1, 2, 3, 4, 5), ("bbb", (0, 1, 2, 3)))
    actual_result = unpack_generators(valid_test_data)
    assert expected_result == actual_result

PROT_COMPARE_DATA = pd.DataFrame.from_dict({
    "2": {"loop_1": "lspsrpgmqd", "loop_2": "fqfhsqymkr"},  # Lowercase
    "5": {"loop_1": "IYCSFVEM", "loop_2": "PQIMNHIGNQKTREW"},  # Different length loops
    "13": {"loop_1": "HWYSFNKKWK", "loop_2": "TVHMNPNKWA", "loop_3": "VEELGPWITV"},  # Three loops
    "142": {"loop_1": "PTQVSEFTRC"}},  # Only one loop
    orient="index")

def test_full_prot_comparison():
    expected_result_1 = ("test_data ['loop_1'] shift(['loop_1'], 'first', 3)",
                         ('2', '5', '13', '142'),
                         (
                          ("2", (81.0, 43.1, 54.35, 72.95)),
                          ("5", (54.35, 66.55, 53.3, 57.8)),
                          ("13", (50.55, 47.4, 78.4, 59.3)),
                          ("142", (64.75, 58.7, 60.7, 77.2))
                         ))
    actual_result_1 = unpack_generators(full_prot_comparison(("test_data", PROT_COMPARE_DATA),
                                                              shift=(["loop_1"], "first", 3)))
    assert expected_result_1 == actual_result_1

    expected_result_2 = ("test_data ['loop_1'] shift(['loop_1'], 'first', 3)",
                         ('2', '5', '13', '142'),
                         (
                          ("13", (50.55, 47.4, 78.4, 59.3)),
                         ))
    actual_result_2 = unpack_generators(full_prot_comparison(("test_data", PROT_COMPARE_DATA),
                                                                shift=(["loop_1"], "first", 3),
                                                                proteins_of_interest=["13"]))
    assert expected_result_2 == actual_result_2

def test_one_prot_all_shifts():
    expected_result = ("test_data protein 142 all shifts", ('2', '5', '13', '142'),
                       (
                        (("('loop_1', 'first', 0)", (64.45, 61.9, 63.0, 56.7)),
                        ("('loop_1', 'first', 1)", (64.75, 61.9, 61.5, 63.6)),
                        ("('loop_1', 'first', 2)", (66.75, 58.7, 61.5, 68.4)),
                        ("('loop_1', 'first', 3)", (64.75, 58.7, 60.7, 77.2)),
                        ("('loop_1', 'first', 4)", (65.95, 56.9, 58.9, 83.8)),
                        ("('loop_1', 'first', 5)", (66.75, 55.7, 55.7, 91.0)),
                        ("('loop_1', 'first', 6)", (66.75, 44.9, 52.5, 100.6)),
                        ("('loop_1', 'first', 7)", (61.05, 40.4, 52.5, 111.4)),
                        ("('loop_1', 'second', 0)", (74.75, 58.5, 64.3, 56.7)),
                        ("('loop_1', 'second', 1)", (74.45, 58.2, 62.3, 63.6)),
                        ("('loop_1', 'second', 2)", (75.65, 57.0, 61.1, 68.4)),
                        ("('loop_1', 'second', 3)", (72.95, 57.8, 59.3, 77.2)),
                        ("('loop_1', 'second', 4)", (64.95, 49.8, 51.3, 83.8)),
                        ("('loop_1', 'second', 5)", (59.75, 48.6, 50.1, 91.0)),
                        ("('loop_1', 'second', 6)", (58.25, 45.6, 54.1, 100.6)),
                        ("('loop_1', 'second', 7)", (57.05, 44.8, 53.3, 111.4))
                       )))
    actual_result = unpack_generators(one_prot_all_shifts("142", ("test_data", PROT_COMPARE_DATA)))
    assert expected_result == actual_result

def test_all_prots_all_shifts():
    test_data = pd.DataFrame.from_dict({
        "prot_1": {"loop_1": "IYCSFVEM", "loop_2": "HIGNQYWMEW"},
        "prot_2": {"loop_1": "MLLSVPLLGAVAE", "loop_2": "WDKPEHIPDPDA"},
        "prot_3": {"loop_1": "PDPSIYAYDNF", "loop_2": "AYAEEFGNETWGVT"}},
        orient="index")
    expected_result = ('highest scores from all_prots_all_shifts, test_data',
                       ('prot_1', 'prot_2', 'prot_3'),
                       (
                        ('prot_1', (234.75, 135.4, 140.55)),
                        ('prot_2', (135.4, 299.85, 164.25)),
                        ('prot_3', (140.55, 164.25, 306.6))
                       ))
    actual_result = unpack_generators(all_prots_all_shifts(("test_data", test_data)))
    assert expected_result == actual_result

def test_amino_acid_properties_matrix():
    # pylint:disable=line-too-long
    scores = (
        ('A', (13.2, 7.85, 4.65, 3.45, 4.95, 8.85, 4.45, 5.7, 3.45, 5.7, 5.7, 5.65, 4.45, 8.4, 2.25, 7.85, 6.65, 7.65, 3.75, 4.5)),
        ('C', (7.85, 15.0, 6.4, 5.2, 3.2, 8.6, 6.2, 3.95, 5.2, 3.95, 3.95, 7.4, 6.2, 6.65, 4.0, 9.6, 8.4, 5.9, 2.0, 2.75)),
        ('D', (4.65, 6.4, 13.2, 10.2, 2.4, 5.4, 7.4, 3.15, 8.4, 3.15, 3.15, 11.3, 7.4, 5.85, 7.2, 6.4, 7.6, 5.1, 1.2, 1.95)),
        ('E', (3.45, 5.2, 10.2, 13.2, 3.6, 4.2, 8.6, 4.35, 9.6, 4.35, 4.35, 7.4, 11.3, 4.65, 8.4, 5.2, 6.4, 3.9, 2.4, 3.15)),
        ('F', (4.95, 3.2, 2.4, 3.6, 13.2, 4.2, 6.4, 8.85, 3.6, 8.85, 8.85, 3.4, 4.6, 6.15, 2.4, 3.2, 4.4, 6.9, 10.2, 10.35)),
        ('G', (8.85, 8.6, 5.4, 4.2, 4.2, 15.0, 5.2, 4.95, 4.2, 4.95, 4.95, 6.4, 5.2, 7.65, 3.0, 8.6, 7.4, 6.9, 3.0, 3.75)),
        ('H', (4.45, 6.2, 7.4, 8.6, 6.4, 5.2, 13.2, 5.35, 10.4, 5.35, 5.35, 8.4, 9.6, 5.65, 9.2, 6.2, 7.4, 4.9, 5.2, 5.95)),
        ('I', (5.7, 3.95, 3.15, 4.35, 8.85, 4.95, 5.35, 13.2, 4.35, 11.4, 9.6, 4.15, 5.35, 6.9, 3.15, 3.95, 6.95, 9.45, 7.65, 8.4)),
        ('K', (3.45, 5.2, 8.4, 9.6, 3.6, 4.2, 10.4, 4.35, 13.2, 4.35, 4.35, 7.4, 8.6, 4.65, 10.2, 5.2, 6.4, 3.9, 2.4, 3.15)),
        ('L', (5.7, 3.95, 3.15, 4.35, 8.85, 4.95, 5.35, 11.4, 4.35, 13.2, 9.6, 4.15, 5.35, 6.9, 3.15, 3.95, 6.95, 9.45, 7.65, 8.4)),
        ('M', (5.7, 3.95, 3.15, 4.35, 8.85, 4.95, 5.35, 9.6, 4.35, 9.6, 13.2, 4.15, 5.35, 6.9, 3.15, 3.95, 5.15, 7.65, 7.65, 8.4)),
        ('N', (5.65, 7.4, 11.3, 7.4, 3.4, 6.4, 8.4, 4.15, 7.4, 4.15, 4.15, 13.2, 8.4, 6.85, 6.2, 7.4, 8.6, 6.1, 2.2, 2.95)),
        ('Q', (4.45, 6.2, 7.4, 11.3, 4.6, 5.2, 9.6, 5.35, 8.6, 5.35, 5.35, 8.4, 13.2, 5.65, 7.4, 6.2, 7.4, 4.9, 3.4, 4.15)),
        ('P', (8.4, 6.65, 5.85, 4.65, 6.15, 7.65, 5.65, 6.9, 4.65, 6.9, 6.9, 6.85, 5.65, 15.0, 3.45, 6.65, 7.85, 8.85, 4.95, 5.7)),
        ('R', (2.25, 4.0, 7.2, 8.4, 2.4, 3.0, 9.2, 3.15, 10.2, 3.15, 3.15, 6.2, 7.4, 3.45, 13.2, 4.0, 5.2, 2.7, 3.6, 4.35)),
        ('S', (7.85, 9.6, 6.4, 5.2, 3.2, 8.6, 6.2, 3.95, 5.2, 3.95, 3.95, 7.4, 6.2, 6.65, 4.0, 13.2, 8.4, 5.9, 2.0, 2.75)),
        ('T', (6.65, 8.4, 7.6, 6.4, 4.4, 7.4, 7.4, 6.95, 6.4, 6.95, 5.15, 8.6, 7.4, 7.85, 5.2, 8.4, 13.2, 8.9, 3.2, 3.95)),
        ('V', (7.65, 5.9, 5.1, 3.9, 6.9, 6.9, 4.9, 9.45, 3.9, 9.45, 7.65, 6.1, 4.9, 8.85, 2.7, 5.9, 8.9, 13.2, 5.7, 6.45)),
        ('W', (3.75, 2.0, 1.2, 2.4, 10.2, 3.0, 5.2, 7.65, 2.4, 7.65, 7.65, 2.2, 3.4, 4.95, 3.6, 2.0, 3.2, 5.7, 13.2, 10.65)),
        ('Y', (4.5, 2.75, 1.95, 3.15, 10.35, 3.75, 5.95, 8.4, 3.15, 8.4, 8.4, 2.95, 4.15, 5.7, 4.35, 2.75, 3.95, 6.45, 10.65, 13.2)))
    expected_result = ("amino acid comparison matrix", ("A", "C", "D", "E", "F", "G", "H", "I", "K",
                       "L", "M", "N", "Q", "P", "R", "S", "T", "V", "W", "Y"), scores)
    actual_result = unpack_generators(amino_acid_properties_matrix())
    assert expected_result == actual_result

def test_compare_peptides():
    expected_result_1 = 58.65
    actual_result_1 = compare_peptides("LSPSRPGMQD", "PTQVSEFTRC", "", "")
    assert expected_result_1 == actual_result_1

    expected_result_2 = 64.75
    actual_result_2 = compare_peptides("LSPSRPGMQD", "PTQVSEFTRC", "second", 3)
    assert expected_result_2 == actual_result_2

def test_find_structural_motifs():
    test_data = (("ILLKCLLWHLEHAVDCTDGTSKPKT", ["Multiple beta-branched residues"]),
                 ("QYHAPWHPKSHKCQQCSYNHAGWVA", ["Big Box N-cap", "Broad Shellman C-cap"]),
                 ("EAAERECYNPKEWFMTACEPWPHFG", ["Box N-cap with a proline", "i-->i+4 Salt bridge"]),
                 ("AHDQKESEIGNWGNFPHFHFSMKGR", ['i-->i+3 Salt bridge', 'AlphaL C-cap']))
    for peptide, expected_result in test_data:
        actual_result = find_structural_motifs(peptide)
        assert expected_result == actual_result
