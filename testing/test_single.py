from pyguide import guide
import os
import glob


def test_single_i():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    file_path_mouse = os.path.join(file_path_1, "example", "gene_list_mouse.txt")
    gene_list = guide.read_gene_list(file_path_2)
    gene_list_mouse = guide.read_gene_list(file_path_mouse)
    assert gene_list[0] == "STAT3"
    guide.order_guides([],
                       gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       check_db=True)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       check_db=False)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list_mouse,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       check_db=True)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list_mouse,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       check_db=False)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)


def test_single_a():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    file_path_mouse = os.path.join(file_path_1, "example", "gene_list_mouse.txt")
    gene_list = guide.read_gene_list(file_path_2)
    gene_list_mouse = guide.read_gene_list(file_path_mouse)
    assert gene_list[0] == "STAT3"
    guide.order_guides([],
                       gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       check_db=True)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       check_db=False)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list_mouse,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       check_db=True)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
    guide.order_guides([],
                       gene_list_mouse,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="single",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       check_db=False)
    single_files = glob.glob(os.path.join(file_path_1, "example", "order_single_Test_*.csv"))
    assert len(single_files) >= 1, "order_guides(single) produced no output CSV"
    for f in single_files + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
        os.remove(f)
