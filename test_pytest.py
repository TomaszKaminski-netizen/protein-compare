"""This module is used for unit testing the protein_compare module."""

from pytest import raises
from protein_compare import one_prot_all_shifts, unpack_generators, full_prot_comparison, \
    compare_peptides, amino_acid_properties_matrix, find_structural_motifs
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

PROT_COMPARE_DATA = {
    "2": {"loop_1": "lspsrpgmqd", "loop_2": "fqfhsqymkr"},  # Lowercase
    "5": {"loop_1": "IYCSFVEM", "loop_2": "PQIMNHIGNQKTREW"},  # Different length loops
    "13": {"loop_1": "HWYSFNKKWK", "loop_2": "TVHMNPNKWA", "loop_3": "VEELGPWITV"},  # Three loops
    "142": {"loop_1": "PTQVSEFTRC"}}  # Only one loop

def test_full_prot_comparison():
    expected_result_1 = ("test_data ['loop_1'] shift(['loop_1'], 'first', 3)",
                         ('2', '5', '13', '142'),
                         (
                          (('81', '43', '54', '73'), '2'),
                          (('54', '67', '53', '58'), '5'),
                          (('51', '47', '78', '59'), '13'),
                          (('65', '59', '61', '77'), '142')
                         ))
    actual_result_1 = unpack_generators(full_prot_comparison(("test_data", PROT_COMPARE_DATA),
                                                              shift=(["loop_1"], "first", 3)))
    assert expected_result_1 == actual_result_1

    expected_result_2 = ("test_data ['loop_1'] shift(['loop_1'], 'first', 3)",
                         ('2', '5', '13', '142'),
                         (
                          (('51', '47', '78', '59'), '13'),
                         ))
    actual_result_2 = unpack_generators(full_prot_comparison(("test_data", PROT_COMPARE_DATA),
                                                                shift=(["loop_1"], "first", 3),
                                                                proteins_of_interest=["13"]))
    assert expected_result_2 == actual_result_2

def test_one_prot_all_shifts():
    expected_result = ("test_data protein 142 all shifts", ('2', '5', '13', '142'),
                       (
                        (('64', '62', '63', '57'), "('loop_1', 'first', 0)"),
                        (('65', '62', '62', '64'), "('loop_1', 'first', 1)"),
                        (('67', '59', '62', '68'), "('loop_1', 'first', 2)"),
                        (('65', '59', '61', '77'), "('loop_1', 'first', 3)"),
                        (('66', '57', '59', '84'), "('loop_1', 'first', 4)"),
                        (('67', '56', '56', '91'), "('loop_1', 'first', 5)"),
                        (('67', '45', '52', '101'), "('loop_1', 'first', 6)"),
                        (('61', '40', '52', '111'), "('loop_1', 'first', 7)"),
                        (('75', '58', '64', '57'), "('loop_1', 'second', 0)"),
                        (('74', '58', '62', '64'), "('loop_1', 'second', 1)"),
                        (('76', '57', '61', '68'), "('loop_1', 'second', 2)"),
                        (('73', '58', '59', '77'), "('loop_1', 'second', 3)"),
                        (('65', '50', '51', '84'), "('loop_1', 'second', 4)"),
                        (('60', '49', '50', '91'), "('loop_1', 'second', 5)"),
                        (('58', '46', '54', '101'), "('loop_1', 'second', 6)"),
                        (('57', '45', '53', '111'), "('loop_1', 'second', 7)")
                       ))
    actual_result = unpack_generators(one_prot_all_shifts("142", ("test_data", PROT_COMPARE_DATA)))
    assert expected_result == actual_result

def test_amino_acid_properties_matrix():
    # pylint:disable=bad-whitespace
    scores = (
        (('13','8','5','3','5','9','4','6','3','6','6','6','4','8','2','8','7','8','4','4'),'A'),
        (('8','15','6','5','3','9','6','4','5','4','4','7','6','7','4','10','8','6','2','3'),'C'),
        (('5','6','13','10','2','5','7','3','8','3','3','11','7','6','7','6','8','5','1','2'),'D'),
        (('3','5','10','13','4','4','9','4','10','4','4','7','11','5','8','5','6','4','2','3'),'E'),
        (('5','3','2','4','13','4','6','9','4','9','9','3','5','6','2','3','4','7','10','10'),'F'),
        (('9','9','5','4','4','15','5','5','4','5','5','6','5','8','3','9','7','7','3','4'),'G'),
        (('4','6','7','9','6','5','13','5','10','5','5','8','10','6','9','6','7','5','5','6'),'H'),
        (('6','4','3','4','9','5','5','13','4','11','10','4','5','7','3','4','7','9','8','8'),'I'),
        (('3','5','8','10','4','4','10','4','13','4','4','7','9','5','10','5','6','4','2','3'),'K'),
        (('6','4','3','4','9','5','5','11','4','13','10','4','5','7','3','4','7','9','8','8'),'L'),
        (('6','4','3','4','9','5','5','10','4','10','13','4','5','7','3','4','5','8','8','8'),'M'),
        (('6','7','11','7','3','6','8','4','7','4','4','13','8','7','6','7','9','6','2','3'),'N'),
        (('4','6','7','11','5','5','10','5','9','5','5','8','13','6','7','6','7','5','3','4'),'Q'),
        (('8','7','6','5','6','8','6','7','5','7','7','7','6','15','3','7','8','9','5','6'),'P'),
        (('2','4','7','8','2','3','9','3','10','3','3','6','7','3','13','4','5','3','4','4'),'R'),
        (('8','10','6','5','3','9','6','4','5','4','4','7','6','7','4','13','8','6','2','3'),'S'),
        (('7','8','8','6','4','7','7','7','6','7','5','9','7','8','5','8','13','9','3','4'),'T'),
        (('8','6','5','4','7','7','5','9','4','9','8','6','5','9','3','6','9','13','6','6'),'V'),
        (('4','2','1','2','10','3','5','8','2','8','8','2','3','5','4','2','3','6','13','11'),'W'),
        (('4','3','2','3','10','4','6','8','3','8','8','3','4','6','4','3','4','6','11','13'),'Y'))
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
