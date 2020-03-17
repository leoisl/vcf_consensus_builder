import unittest
from io import StringIO
from vcf_consensus_builder.vcf_io import read_vcf
from vcf_consensus_builder.vcf_consensus_builder_core import get_interval_tree_for_vcf_of_a_single_chrom, InconsistentVCFException, MultipleChromException
from intervaltree import Interval, IntervalTree


class Test_get_gt_from_sample_info(unittest.TestCase):
    def test_get_interval_tree_for_vcf___vcf_with_duplicated_record_should_fail(self):
        vcf = StringIO(
"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	2	.	G	T	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
""")
        df = read_vcf(vcf)
        self.assertRaises(InconsistentVCFException, get_interval_tree_for_vcf_of_a_single_chrom, df)


    def test_get_interval_tree_for_vcf___should_be_ok(self):
        vcf = StringIO(
"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	3	.	ACCGT	CCCC	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	10	.	GGA	TTT	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
""")
        df = read_vcf(vcf)
        actual = get_interval_tree_for_vcf_of_a_single_chrom(df)
        expected = IntervalTree([Interval(2, 3, ('G', 'A')), Interval(3, 8, ('ACCGT', 'CCCC')), Interval(10, 13, ('GGA', 'TTT'))])
        self.assertEqual(actual, expected)


    def test_get_interval_tree_for_vcf___two_chrom___should_fail(self):
        vcf = StringIO(
"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_2	3	.	A	C	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
""")
        df = read_vcf(vcf)
        self.assertRaises(MultipleChromException, get_interval_tree_for_vcf_of_a_single_chrom, df)



    def test_get_interval_tree_for_vcf___contain_duplicates_but_GT_is_null_or_ref_and_are_disconsidered___should_be_ok(self):
        vcf = StringIO(
"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	0/1:255,186,0
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	./1:255,186,0
""")
        df = read_vcf(vcf)
        actual = get_interval_tree_for_vcf_of_a_single_chrom(df)
        expected = IntervalTree([Interval(2, 3, ('G', 'A'))])
        self.assertEqual(actual, expected)


    def test_get_interval_tree_for_vcf___contain_enveloped___should_be_ok(self):
        vcf = StringIO(
"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
ref_1	2	.	G	A	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	3	.	A	C	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	3	.	ACCGT	CCCC	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	3	.	AC	CC	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
ref_1	10	.	GGA	TTT	228	.	DP=9;VDB=0.000651286;SGB=-0.693147;RPB=1;MQB=1;MQSB=0.110884;BQB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,4;MQ=50	GT:PL	1/1:255,186,0
""")
        df = read_vcf(vcf)
        actual = get_interval_tree_for_vcf_of_a_single_chrom(df)
        expected = IntervalTree([Interval(2, 3, ('G', 'A')), Interval(3, 4, ('A', 'C')), Interval(3, 8, ('ACCGT', 'CCCC')),
                                 Interval(3, 5, ('AC', 'CC')), Interval(10, 13, ('GGA', 'TTT'))])
        self.assertEqual(actual, expected)
