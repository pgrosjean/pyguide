from pyguide import guide
import os


def test_arrayed_human_i():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    gene_list = guide.read_gene_list(file_path_2)
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       kampmann_lab=True)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       kampmann_lab=False)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       kampmann_lab=True)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       kampmann_lab=False)


def test_arrayed_order_a():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    gene_list = guide.read_gene_list(file_path_2)
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       kampmann_lab=True)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="human",
                       kampmann_lab=False)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       kampmann_lab=True)
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="arrayed",
                       base_dir=os.path.join(file_path_1, "example"),
                       organism="mouse",
                       kampmann_lab=False)
