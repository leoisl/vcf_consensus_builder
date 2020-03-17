import unittest
from io import StringIO
from vcf_consensus_builder.vcf_consensus_builder_core import ensure_there_are_no_overlapping_records, InconsistentVCFException
from intervaltree import Interval, IntervalTree


class Test_ensure_there_are_no_overlapping_records(unittest.TestCase):
    def test_ensure_there_are_no_overlapping_records___no_overlapping_records___far_away(self):
        interval_tree = IntervalTree([Interval(1, 5), Interval(10, 15), Interval(20, 25)])
        ensure_there_are_no_overlapping_records(interval_tree)
        self.assertTrue(True)

    def test_ensure_there_are_no_overlapping_records___no_overlapping_records___close(self):
        interval_tree = IntervalTree([Interval(1, 5), Interval(5, 10), Interval(10, 15)])
        ensure_there_are_no_overlapping_records(interval_tree)
        self.assertTrue(True)

    def test_ensure_there_are_no_overlapping_records___one_overlap(self):
        interval_tree = IntervalTree([Interval(1, 3), Interval(5, 10), Interval(9, 15)])
        self.assertRaises(InconsistentVCFException, ensure_there_are_no_overlapping_records, interval_tree)

    def test_ensure_there_are_no_overlapping_records___one_overlap_equals(self):
        interval_tree = IntervalTree([Interval(1, 3), Interval(5, 10, 0), Interval(5, 10, 1)])
        self.assertRaises(InconsistentVCFException, ensure_there_are_no_overlapping_records, interval_tree)

    def test_ensure_there_are_no_overlapping_records___one_overlap_enveloped(self):
        interval_tree = IntervalTree([Interval(1, 3), Interval(5, 8, 0), Interval(5, 10, 1)])
        self.assertRaises(InconsistentVCFException, ensure_there_are_no_overlapping_records, interval_tree)
