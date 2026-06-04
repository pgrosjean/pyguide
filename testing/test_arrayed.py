from pyguide import guide
import os
import glob


def test_arrayed_human_i():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    file_path_mouse = os.path.join(file_path_1, "example", "gene_list_mouse.txt")
    gene_list = guide.read_gene_list(file_path_2)
    gene_list_mouse = guide.read_gene_list(file_path_mouse)
    assert gene_list[0] == "STAT3"
    try:
        guide.order_guides([],
                           gene_list,
                           name="Test",
                           ai_status="i",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="human",
                           check_db=True)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list,
                           name="Test",
                           ai_status="i",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="human",
                           check_db=False)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list_mouse,
                           name="Test",
                           ai_status="i",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="mouse",
                           check_db=True)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list_mouse,
                           name="Test",
                           ai_status="i",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="mouse",
                           check_db=False)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)


def test_arrayed_order_a():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    file_path_mouse = os.path.join(file_path_1, "example", "gene_list_mouse.txt")
    gene_list = guide.read_gene_list(file_path_2)
    gene_list_mouse = guide.read_gene_list(file_path_mouse)
    assert gene_list[0] == "STAT3"
    try:
        guide.order_guides([],
                           gene_list,
                           name="Test",
                           ai_status="a",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="human",
                           check_db=True)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list,
                           name="Test",
                           ai_status="a",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="human",
                           check_db=False)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list_mouse,
                           name="Test",
                           ai_status="a",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="mouse",
                           check_db=True)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
    try:
        guide.order_guides([],
                           gene_list_mouse,
                           name="Test",
                           ai_status="a",
                           guides_per_gene=5,
                           order_format="arrayed",
                           base_dir=os.path.join(file_path_1, "example"),
                           organism="mouse",
                           check_db=False)
        arrayed_files = glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv"))
        assert len(arrayed_files) >= 1, "order_guides(arrayed) produced no output CSV"
    finally:
        for f in glob.glob(os.path.join(file_path_1, "example", "order_arrayed_Test_*.csv")) \
                + glob.glob(os.path.join(file_path_1, "example", "log_file_*Test_*.txt")):
            if os.path.exists(f):
                os.remove(f)
