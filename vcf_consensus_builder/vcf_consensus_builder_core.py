"""Main module."""
import logging
from io import TextIOWrapper

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from intervaltree import Interval, IntervalTree

from vcf_consensus_builder.vcf_io import read_vcf

from typing import List

import re

logger = logging.getLogger(__name__)


def mutate_seqs_according_given_positions(seq_records: List[SeqRecord], mutable_seqs: List[MutableSeq],
                                          contig_and_position: pd.DataFrame, replacement_char: str) -> None:
    for seq_record, mutable_seq in zip(seq_records, mutable_seqs):
        positions_for_this_seq = contig_and_position[contig_and_position.contig == seq_record.id].position
        for position in positions_for_this_seq:
            mutable_seq[position] = replacement_char



def replace_low_depth_positions(ref_fasta: str,
                                depths_file: str,
                                no_coverage: int = 0,
                                low_coverage: int = 5,
                                no_cov_char: str = '-',
                                low_cov_char: str = 'N',
                                *args,
                                **kwargs) -> List[SeqRecord]:
    """Replace no and low coverage depth positions with specified characters.

    Args:
        ref_fasta: FASTA file path
        depths_file: Samtools depth tab-delimited file path
        no_coverage: No coverage depth threshold. At and below this coverage positions are changed to the `no_cov_char`
        low_coverage: Low coverage depth threshold. Below this coverage positions are changed to the `low_cov_char`
        no_cov_char: Character to substitute at no coverage positions (default: "-")
        low_cov_char: Character to substitute at low coverage positions (default: "N")

    Returns:
        SeqRecord with low coverage depth positions substituted with `unmapped_char`
    """
    assert len(no_cov_char) == 1, '"no_cov_char" must be a str of length of 1, e.g. "-"/single dash character.'
    assert len(low_cov_char) == 1, '"low_cov_char" must be a str of length of 1, e.g. "-"/single dash character.'
    assert no_coverage <= low_coverage, '"low_coverage" must be greater than or equal to "no_coverage"'
    seq_records: List[SeqRecord] = list(SeqIO.parse(ref_fasta,
                         format='fasta'))
    logger.debug(f'Number of contigs: {len(seq_records)}')
    df: pd.DataFrame = pd.read_csv(depths_file,
                                   header=None,
                                   names=['contig', 'position', 'coverage'], sep='\t')
    df.position = df.position - 1 # 0-based indexing fix

    no_coverage_positions: pd.DataFrame = df[df.coverage <= no_coverage][["contig", "position"]]
    logger.info(f'No ({no_coverage}X) coverage positions: {len(no_coverage_positions)}')
    low_coverage_positions: pd.DataFrame = df[(df.coverage > no_coverage) & (df.coverage < low_coverage)][["contig", "position"]]
    logger.info(f'Low (<{low_coverage}X) coverage positions: {len(low_coverage_positions)}')
    if len(low_coverage_positions) == 0 and len(no_coverage_positions) == 0:
        logger.info(f'No positions with low (<{low_coverage}X) or no ({no_coverage}X) coverage. '
                    f'No need to mask any positions in the reference sequence')
        return seq_records
    mutable_seqs: List[MutableSeq] = [seq_record.seq.tomutable() for seq_record in seq_records]
    mutate_seqs_according_given_positions(seq_records, mutable_seqs, low_coverage_positions, low_cov_char)
    mutate_seqs_according_given_positions(seq_records, mutable_seqs, no_coverage_positions, no_cov_char)

    mutated_seq_records = [SeqRecord(seq=mutable_seq.toseq(), id=seq_record.id, name=seq_record.name, description=seq_record.description)
                           for mutable_seq, seq_record in zip(mutable_seqs, seq_records)]
    return mutated_seq_records


def consensus_segment(seq: str,
                      curr_position: int,
                      ref_variant: str,
                      alt_variant: str,
                      prev_position: int = 0) -> (str, int):
    segment_before_curr_position = seq[prev_position:(curr_position - 1)]
    segment_on_curr_position = seq[(curr_position - 1):]
    ref_variant_is_correct = segment_on_curr_position.startswith(ref_variant)
    if not ref_variant_is_correct:
        raise Exception(f"Error: {ref_variant} expected in the start of {segment_on_curr_position}")

    segment_seq = segment_before_curr_position + alt_variant
    next_position = curr_position + len(ref_variant) - 1
    return segment_seq, next_position


def get_gt_from_sample_info(sample_info: str):
    gt_as_string = re.search(r'^(.+?)[/:]??', sample_info).group(1)
    try:
        gt_as_int = int(gt_as_string)
        return gt_as_int
    except ValueError:
        return -1

class InconsistentVCFException(Exception):
    pass

class MultipleChromException(Exception):
    pass

def get_interval_tree_for_vcf_of_a_single_chrom(df_vcf):
    only_a_single_chrom_in_the_df = len(df_vcf["CHROM"].unique()) == 1
    if not only_a_single_chrom_in_the_df:
        raise MultipleChromException


    interval_tree = IntervalTree()
    intervals_already_inserted = set() # this is done because it is easier to query than interval_tree for exact intervals and without data
    for _, curr_var in df_vcf.iterrows():
        start_pos = curr_var.POS
        ref = curr_var.REF
        end_pos = start_pos + len(ref)
        sample_info = curr_var[-1]
        GT = get_gt_from_sample_info(sample_info)
        variant_should_not_be_applied = GT==-1 or GT==0
        if variant_should_not_be_applied:
            continue

        alleles = [curr_var.REF] + curr_var.ALT.split(",")
        alt = alleles[GT]
        data = (ref, alt)
        interval_already_in_tree = (start_pos, end_pos) in intervals_already_inserted
        if interval_already_in_tree:
            raise InconsistentVCFException(f"[FATAL] There are more than one VCF record with POS = {start_pos} and REF = {ref}")
        interval_tree[start_pos:end_pos] = data
        intervals_already_inserted.add((start_pos, end_pos))

    return interval_tree


def ensure_there_are_no_overlapping_records(tree):
    for interval in tree:
        overlaps = tree.overlap(interval.begin, interval.end)
        does_not_overlap_with_any_other_interval = len(overlaps) == 1
        if not does_not_overlap_with_any_other_interval:
            raise InconsistentVCFException(f"{interval} overlaps with {overlaps}")


def get_super_records_from_interval_tree(interval_tree):
    all_envelopped_intervals = IntervalTree()
    for interval in interval_tree:
        for enveloped_interval in interval_tree.envelop(interval):
            if interval != enveloped_interval:
                envelopped_interval_is_consistent = \
                    interval.begin == enveloped_interval.begin and \
                    interval.end > enveloped_interval.end and \
                    interval.data[0].startswith(enveloped_interval.data[0]) and \
                    interval.data[1].startswith(enveloped_interval.data[1])

                if not envelopped_interval_is_consistent:
                    raise InconsistentVCFException(f"[FATAL]: {enveloped_interval} is not consistent with {interval}")
                all_envelopped_intervals.add(enveloped_interval)

    super_records = interval_tree - all_envelopped_intervals
    return super_records


def get_records_to_be_applied(df_vcf):
    interval_tree_with_valid_records = get_interval_tree_for_vcf_of_a_single_chrom(df_vcf)
    interval_tree_with_records_to_be_applied = get_super_records_from_interval_tree(interval_tree_with_valid_records)
    ensure_there_are_no_overlapping_records(interval_tree_with_records_to_be_applied)
    records_to_be_applied = []
    for interval in sorted(interval_tree_with_records_to_be_applied):
        records_to_be_applied.append({
            "begin": interval.begin,
            "end": interval.end,
            "ref": interval.data[0],
            "alt": interval.data[1]
        })
    return records_to_be_applied


def create_consensus_sequences(ref_seq_records: List[SeqRecord], df_vcf_tsv: pd.DataFrame) -> List[str]:
    consensus_sequences: List[str] = []
    for ref_seq_record in ref_seq_records:
        df_vcf_tsv_for_this_record = df_vcf_tsv[df_vcf_tsv.CHROM == ref_seq_record.id]
        records_to_be_applied = get_records_to_be_applied(df_vcf_tsv_for_this_record)
        consensus_sequence = create_cons_seq(str(ref_seq_record.seq), records_to_be_applied)
        consensus_sequences.append(consensus_sequence)
    return consensus_sequences

def create_cons_seq(seq: str, records_to_be_applied: List) -> str:
    """Create consensus sequence given a reference sequence and a table of a Snippy vcf_to_tab

    Args:
        seq: Reference sequence
        df_vcf: VCF file DataFrame

    Returns:
        Consensus sequence
    """
    segments = []
    prev_position = 0
    for record in records_to_be_applied:
        segment, prev_position = consensus_segment(seq=seq,
                                                   curr_position=record["begin"],
                                                   ref_variant=record["ref"],
                                                   alt_variant=record["alt"],
                                                   prev_position=prev_position)
        segments.append(segment)
    # append the rest of the reference sequence
    segments.append(seq[prev_position:])
    return ''.join(segments)


def consensus(ref_fasta,
              vcf_file,
              depths_file,
              output_fasta,
              no_coverage: int = 0,
              low_coverage: int = 5,
              no_cov_char: str = '-',
              low_cov_char: str = 'N',
              *args,
              **kwargs):
    ref_seq_records: List[SeqRecord] = replace_low_depth_positions(**locals())
    with open(vcf_file) as fh:
        df_vcf_tsv: pd.DataFrame = read_vcf(fh)
    logger.debug(f'df_vcf_tsv shape: {df_vcf_tsv.shape}')

    consensus_seqs: List[str] = create_consensus_sequences(ref_seq_records, df_vcf_tsv)
    with open(output_fasta, 'w') if not isinstance(output_fasta, TextIOWrapper) else output_fasta as f:
        for ref_seq_record, consensus_seq in zip(ref_seq_records, consensus_seqs):
            print(f'>{ref_seq_record.description}', file=f)
            print(consensus_seq, file=f)
