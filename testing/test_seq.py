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


def test_write_pooled_seq_log_file(tmp_path):
    primer_df = pd.DataFrame({
        'guide_id': ['G1', 'G2'],
        'seq': ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
        'left_primers': ['LEFTSEQ', 'LEFTSEQ'],
        'right_primers': ['RIGHTSEQ', 'RIGHTSEQ'],
        'lib_num': [0, 0],
    })
    guide.write_pooled_seq_log_file("TestUser", str(tmp_path), primer_df)
    log_files = list(tmp_path.glob("log_file_pooled_seq_TestUser_*.txt"))
    assert len(log_files) == 1
    content = log_files[0].read_text()
    assert "TestUser" in content
    assert "LEFTSEQ" in content
    assert "RIGHTSEQ" in content


def test_collate_seq_and_order():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    seq_file = os.path.join(file_path_1, "example", "seq_guide_list.txt")
    example_dir = os.path.join(file_path_1, "example")

    try:
        # Step 1: collate
        collate_seq.main(["--sequence_files", seq_file])
        collated = glob.glob(os.path.join(example_dir, "collated_seq_wishlist_*.txt"))
        assert len(collated) >= 1, "collate_seq produced no collated file"
        collated_file = sorted(collated, key=os.path.getmtime)[-1]

        # Step 2: read collated file and build primer_df
        names, seqs, left_primers, right_primers, lib_nums = guide.read_gene_list_pooled_seq(collated_file)
        primer_df = pd.DataFrame({
            'guide_id': names,
            'seq': seqs,
            'left_primers': left_primers,
            'right_primers': right_primers,
            'lib_num': lib_nums,
        })

        # Step 3: order
        guide.order_guides(
            guide_ids=[],
            gene_names=[],
            name="Test",
            ai_status="i",
            guides_per_gene=5,
            order_format="pooled-seq",
            base_dir=example_dir,
            check_db=False,
            organism="human",
            primer_df=primer_df,
        )

        # Step 4: verify output file exists and has content
        order_files = glob.glob(os.path.join(example_dir, "order_pooled_Test_*.txt"))
        assert len(order_files) >= 1, "order_guides(pooled-seq) produced no output TXT"
        with open(sorted(order_files, key=os.path.getmtime)[-1]) as fh:
            lines = [l for l in fh.readlines() if l.strip()]
        assert len(lines) >= 6, f"Expected at least 6 lines (2 per guide x 3 guides), got {len(lines)}"

    finally:
        for f in (
            glob.glob(os.path.join(example_dir, "order_pooled_Test_*.txt"))
            + glob.glob(os.path.join(example_dir, "order_pooled_Test_*_info.csv"))
            + glob.glob(os.path.join(example_dir, "log_file_pooled_seq_Test_*.txt"))
            + glob.glob(os.path.join(example_dir, "collated_seq_wishlist_*.txt"))
        ):
            if os.path.exists(f):
                os.remove(f)
