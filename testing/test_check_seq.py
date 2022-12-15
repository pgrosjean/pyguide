from pyguide import check_seq


def check_collate_files():
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, 'example', 'seq_files')
    collated_files = check_seq.collate_files(file_path, ".seq")
    assert len(collated_files) == 26


def test_check_seq():
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, 'example', 'seq_files')
    collated_files = check_seq.collate_files(file_path, ".seq")
    guide_ids = sorted(check_seq.check_seq(collated_files, "human", "i"))
    assert guide_ids[0] == 'ASB1_-_239335717.23-P1P2'
    assert guide_ids[-1] == 'SOCS6_+_67956710.23-P1P2'
