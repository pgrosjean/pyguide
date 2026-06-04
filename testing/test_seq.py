import pytest
import os
import glob
import pandas as pd
from pyguide import collate_seq, guide


def test_validate_sequences_valid(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("MY_GUIDE_1\tACGTACGTACGTACGTACGT\nMY_GUIDE_2\tTGCATGCATGCATGCATGCA\n")
    df = collate_seq.validate_sequence_file(str(f))
    assert len(df) == 2
    assert list(df['name']) == ['MY_GUIDE_1', 'MY_GUIDE_2']
    assert list(df['seq']) == ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA']


def test_validate_sequences_wrong_length(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("SHORT_GUIDE\tACGT\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))


def test_validate_sequences_invalid_chars(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("BAD_CHARS\tACGTACGNACGTACGTACGT\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))


def test_validate_sequences_duplicate_names(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("GUIDE_1\tACGTACGTACGTACGTACGT\nGUIDE_1\tTGCATGCATGCATGCATGCA\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))


def test_read_gene_list_pooled_seq(tmp_path):
    f = tmp_path / "collated.txt"
    f.write_text(
        "GUIDE_A\tACGTACGTACGTACGTACGT\tLEFT1\tRIGHT1\t0\n"
        "GUIDE_B\tTGCATGCATGCATGCATGCA\tLEFT1\tRIGHT1\t0\n"
    )
    names, seqs, lefts, rights, lib_nums = guide.read_gene_list_pooled_seq(str(f))
    assert names == ['GUIDE_A', 'GUIDE_B']
    assert seqs == ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA']
    assert lefts == ['LEFT1', 'LEFT1']
    assert rights == ['RIGHT1', 'RIGHT1']
    assert lib_nums == [0, 0]
