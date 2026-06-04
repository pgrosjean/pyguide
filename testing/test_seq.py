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
