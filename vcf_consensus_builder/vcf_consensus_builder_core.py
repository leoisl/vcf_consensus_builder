"""Main module."""
import logging
from io import TextIOWrapper

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq

from vcf_consensus_builder.vcf_io import (read_vcf, VCF_COL_DTYPES)

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
    segment_seq = seq[prev_position:(curr_position - 1)] + alt_variant
    next_position = curr_position + len(ref_variant) - 1
    if len(alt_variant) != len(ref_variant):
        logger.info(f'At position {curr_position}, ALT="{alt_variant}" '
                    f'(n={len(alt_variant)} VS REF="{ref_variant}" (n={len(ref_variant)})')
        logger.info(f'Previous position={prev_position} | curr_position={curr_position} | next={next_position}')
    return segment_seq, next_position


def create_consensus_sequences(ref_seq_records: List[SeqRecord], df_vcf_tsv: pd.DataFrame) -> List[str]:
    consensus_sequences: List[str] = []
    for ref_seq_record in ref_seq_records:
        df_vcf_tsv_for_this_record = df_vcf_tsv[df_vcf_tsv.CHROM == ref_seq_record.id]
        consensus_sequence = create_cons_seq(str(ref_seq_record.seq), df_vcf_tsv_for_this_record)
        consensus_sequences.append(consensus_sequence)
    return  consensus_sequences


def create_cons_seq(seq: str, df_vcf: pd.DataFrame) -> str:
    """Create consensus sequence given a reference sequence and a table of a Snippy vcf_to_tab

    Args:
        seq: Reference sequence
        df_vcf: VCF file DataFrame

    Returns:
        Consensus sequence
    """
    segments = []
    prev_position = 0
    for _, curr_var in df_vcf.iterrows():
        if prev_position > curr_var.POS - 1:
            logger.warning(
                f'Skipping variant (ALT={curr_var.ALT}) at {curr_var.POS} (previous position ({prev_position}) >= POS)')
            continue

        sample_info = curr_var[-1]
        GT = int(re.search(r'\d+', sample_info).group())
        alleles = [curr_var.REF] + curr_var.ALT.split(",")
        alt = alleles[GT]
        segment, prev_position = consensus_segment(seq=seq,
                                                   curr_position=curr_var.POS,
                                                   ref_variant=curr_var.REF,
                                                   alt_variant=alt,
                                                   prev_position=prev_position)
        segments.append(segment)
    # append the rest of the reference sequence
    segments.append(seq[prev_position:])
    return ''.join(segments)


def consensus(ref_fasta,
              vcf_file,
              depths_file,
              sample_name,
              output_fasta,
              no_coverage: int = 0,
              low_coverage: int = 5,
              no_cov_char: str = '-',
              low_cov_char: str = 'N',
              *args,
              **kwargs):
    ref_seq_records: List[SeqRecord] = replace_low_depth_positions(**locals())
    df_vcf_tsv: pd.DataFrame = read_vcf(vcf_file)
    logger.debug(f'df_vcf_tsv shape: {df_vcf_tsv.shape}')

    consensus_seqs: List[str] = create_consensus_sequences(ref_seq_records, df_vcf_tsv)
    with open(output_fasta, 'w') if not isinstance(output_fasta, TextIOWrapper) else output_fasta as f:
        for ref_seq_record, consensus_seq in zip(ref_seq_records, consensus_seqs):
            print(f'>{ref_seq_record.id}', file=f)
            print(consensus_seq, file=f)
