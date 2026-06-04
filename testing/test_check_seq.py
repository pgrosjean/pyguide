from pyguide import check_seq
import os


def test_collate_files():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "seq_files")
    collated_files = check_seq.collate_files(file_path_2, ".seq")
    assert len(collated_files) == 26


def test_check_seq():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "seq_files")
    collated_files = check_seq.collate_files(file_path_2, ".seq")
    guide_ids = sorted(check_seq.check_seq(collated_files, "human", "i").astype(str))
    assert guide_ids[0] == 'ASB1_-_239335717.23-P1P2'
    assert guide_ids[-1] == 'SOCS6_+_67956710.23-P1P2'


def test_check_seq_crispra():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "seq_files_crispra")
    collated_files = check_seq.collate_files(file_path_2, ".seq")
    # Running with CRISPRa database should not raise
    guide_ids = check_seq.check_seq(collated_files, "human", "a")
    assert guide_ids is not None
    assert len(guide_ids) == len(collated_files)
    assert any(g is not None for g in guide_ids), "CRISPRa check_seq returned all None — lookup may be broken"
