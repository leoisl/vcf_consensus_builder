import unittest
from io import StringIO
from vcf_consensus_builder.vcf_io import read_vcf
from vcf_consensus_builder.vcf_consensus_builder_core import get_super_records_from_interval_tree, InconsistentVCFException
from intervaltree import Interval, IntervalTree


class Test_get_gt_from_sample_info(unittest.TestCase):
    def test_get_super_records_from_interval_tree___no_overlaps___everything_is_a_super_record(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')), Interval(3, 8, ('ACCGT', 'CCCC')), Interval(10, 13, ('GGA', 'TTT'))])
        expected = get_super_records_from_interval_tree(interval_tree)
        actual = interval_tree
        self.assertEqual(actual, expected)

    def test_get_super_records_from_interval_tree___several_overlaps(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(3, 4, ('A', 'C')),
                                      Interval(3, 5, ('AC', 'CC')),
                                      Interval(3, 6, ('ACC', 'CCC')),
                                      Interval(3, 7, ('ACCG', 'CCC')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        expected = IntervalTree([Interval(2, 3, ('G', 'A')), Interval(3, 8, ('ACCGT', 'CCCC')), Interval(10, 13, ('GGA', 'TTT'))])
        actual = get_super_records_from_interval_tree(interval_tree)
        self.assertEqual(actual, expected)

    def test_get_super_records_from_interval_tree___one_overlap_begin_match(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(3, 10, ('ACCGTGG', 'CCCCA')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        expected = IntervalTree(
            [Interval(2, 3, ('G', 'A')), Interval(3, 10, ('ACCGTGG', 'CCCCA')), Interval(10, 13, ('GGA', 'TTT'))])
        actual = get_super_records_from_interval_tree(interval_tree)
        self.assertEqual(actual, expected)


    def test_get_super_records_from_interval_tree___one_overlap_totally_inside(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(4, 7, ('CCG', 'CC')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        expected = IntervalTree(
            [Interval(2, 3, ('G', 'A')), Interval(3, 8, ('ACCGT', 'CCCC')), Interval(10, 13, ('GGA', 'TTT'))])
        actual = get_super_records_from_interval_tree(interval_tree)
        self.assertEqual(actual, expected)

    def test_get_super_records_from_interval_tree___one_overlap_end_match(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(6, 8, ('GT', 'CC')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        expected = IntervalTree(
            [Interval(2, 3, ('G', 'A')), Interval(3, 8, ('ACCGT', 'CCCC')), Interval(10, 13, ('GGA', 'TTT'))])
        actual = get_super_records_from_interval_tree(interval_tree)
        self.assertEqual(actual, expected)

    def test_get_super_records_from_interval_tree___one_overlap_ref_not_consistent_in_the_begin(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(3, 5, ('AG', 'CC')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        self.assertRaises(InconsistentVCFException, get_super_records_from_interval_tree, interval_tree)

    def test_get_super_records_from_interval_tree___one_overlap_ref_not_consistent_in_the_middle(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(5, 6, ('A', 'C')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        self.assertRaises(InconsistentVCFException, get_super_records_from_interval_tree, interval_tree)

    def test_get_super_records_from_interval_tree___one_overlap_ref_not_consistent_in_the_end(self):
        interval_tree = IntervalTree([Interval(2, 3, ('G', 'A')),
                                      Interval(7, 8, ('A', 'C')),
                                      Interval(3, 8, ('ACCGT', 'CCCC')),
                                      Interval(10, 13, ('GGA', 'TTT'))])
        self.assertRaises(InconsistentVCFException, get_super_records_from_interval_tree, interval_tree)
