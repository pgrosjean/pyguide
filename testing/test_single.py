from pyguide import guide
import os


def test_single_order_i():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    gene_list = guide.read_gene_list(file_path_2)
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                 name="Test",
                 ai_status="i",
                 guides_per_gene=5,
                 order_format="single",
                 base_dir=os.path.join(file_path_1, "example"))


def test_single_order_a():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    gene_list = guide.read_gene_list(file_path_2)
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                 name="Test",
                 ai_status="a",
                 guides_per_gene=5,
                 order_format="single",
                 base_dir=os.path.join(file_path_1, "example"))