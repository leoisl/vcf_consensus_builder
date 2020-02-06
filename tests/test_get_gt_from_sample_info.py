import unittest
from vcf_consensus_builder.vcf_consensus_builder_core import get_gt_from_sample_info


class Test_get_gt_from_sample_info(unittest.TestCase):
    def test_correct_gt_with_slash(self):
        sample_info = "3/3:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_without_slash(self):
        sample_info = "3:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_colon_comes_before_slash(self):
        sample_info = "3:rest/asd"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_has_no_colon(self):
        sample_info = "3/3"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_has_no_colon_no_slash(self):
        sample_info = "3"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_several_slashes_and_colons(self):
        sample_info = "3/3:rest/asd:qwe/ewr"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_several_slashes_and_colons_no_slash_in_first(self):
        sample_info = "3:rest/asd:qwe/ewr"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_correct_gt_no_slashes_several_colons(self):
        sample_info = "3:rest:qwe"
        actual = get_gt_from_sample_info(sample_info)
        expected = 3
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_dot_with_slash(self):
        sample_info = "./.:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_dot_without_slash(self):
        sample_info = ".:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_word_with_slash(self):
        sample_info = "asd/qwe:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_word_without_slash(self):
        sample_info = "asd:rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_word_without_colon(self):
        sample_info = "asd/rest"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)

    def test_incorrect_gt_is_word_without_colon_without_slash(self):
        sample_info = "asd"
        actual = get_gt_from_sample_info(sample_info)
        expected = -1
        self.assertEqual(actual, expected)
