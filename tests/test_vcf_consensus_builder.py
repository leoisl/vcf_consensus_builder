#!/usr/bin/env python

"""Tests for `vcf_consensus_builder` package."""

import re
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO
from click.testing import CliRunner

import cli

VCF = 'tests/data/test.vcf'
REF_FASTA = 'tests/data/ref.fa'
DEPTHS = 'tests/data/test-depths.tsv'
OUTPUT_FASTA = 'out.fa'


def test_command_line_interface_help():
    """Test the CLI."""
    runner = CliRunner()
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert re.search(r'--help\s+Show this message and exit.', help_result.output, flags=re.DOTALL) is not None

def test_command_line_interface_vcf_and_ref():
    # Test replacing multi-char deletion and SNPs and setting sample name via command-line
    runner = CliRunner()
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        result = runner.invoke(cli.main,
                               ['-v', VCF,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output])
        assert result.exit_code == 0
        records = list(SeqIO.parse(full_fasta_output, 'fasta'))
        assert str(records[0].seq) == 'NACCGTANACAATAN--'
        assert str(records[1].seq) == 'TATACACATCCACGGC-N-N-N-N-N--T'

def test_command_line_interface_no_coverage_threshold():
    # Test replacing multi-char deletion and SNPs and setting sample name via command-line
    runner = CliRunner()

    # Test changing no coverage threshold
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        result = runner.invoke(cli.main,
                               ['-v', VCF,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output,
                                '--no-coverage', 5])
        assert result.exit_code == 0
        records = list(SeqIO.parse(full_fasta_output, 'fasta'))
        assert str(records[0].seq) == '-ACCGTA-ACAAT----', 'Positions <= 5X coverage must be replaced with "-"'
        assert str(records[1].seq) == 'TATACACATCCACG---------------', 'Positions <= 5X coverage must be replaced with "-"'

def test_command_line_interface_low_and_no_coverage_threshold_with_other_chars():
    # Test replacing low and no coverage characters with other characters than default N and - respectively
    runner = CliRunner()

    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        result = runner.invoke(cli.main,
                               ['-v', VCF,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output,
                                '--no-cov-char', '=',
                                '--low-cov-char', '@'])
        assert result.exit_code == 0
        records = list(SeqIO.parse(full_fasta_output, 'fasta'))
        assert str(records[0].seq) == '@ACCGTA@ACAATA@==', \
            'No coverage positions must be replaced with "=". Low coverage (<5X) positions must be replaced with "@".'
        assert str(records[1].seq) == 'TATACACATCCACGGC=@=@=@=@=@==T', \
            'No coverage positions must be replaced with "=". Low coverage (<5X) positions must be replaced with "@".'
