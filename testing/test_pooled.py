from pyguide import guide, pool
import pandas as pd
import os


def test_pooled_order_i():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    file_path_3 = os.path.join(file_path_1, "example", "gene_list_2.txt")
    file_path_4 = os.path.join(file_path_1, "example", "test_primer_file.txt")
    file_path_5 = os.path.join(file_path_1, "example", "collated_pooled_wish_list_test.txt")
    file_path_mouse = os.path.join(file_path_1, "example", "gene_list_mouse.txt")
    file_path_mouse_2 = os.path.join(file_path_1, "example", "collated_pooled_wish_list_test_mouse.txt")
    # Testing the pool # single and double test
    pool.main(["--wishlist_files", file_path_2, file_path_3])
    pool.main(["--wishlist_files", file_path_2, file_path_3, "--primer_file", file_path_4])
    pool.main(["--wishlist_files", file_path_mouse])

    gene_list, left_primers, right_primers, lib_num = guide.read_gene_list_pooled(file=file_path_5)
    primer_df = pd.DataFrame({'gene_symbol': gene_list,
                              'left_primers': left_primers,
                              'right_primers': right_primers,
                              'lib_num': lib_num})

    glm, lpm, rpm, lnm = guide.read_gene_list_pooled(file=file_path_mouse_2)
    primer_df_mouse = pd.DataFrame({'gene_symbol': glm,
                                    'left_primers': lpm,
                                    'right_primers': rpm,
                                    'lib_num': lnm})
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="human",
                       primer_df=primer_df)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="mouse",
                       primer_df=primer_df_mouse)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=False,
                       organism="mouse",
                       primer_df=primer_df_mouse)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="i",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="mouse",
                       primer_df=primer_df_mouse)


def test_pooled_order_a():
    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_5 = os.path.join(file_path_1, "example", "collated_pooled_wish_list_test.txt")
    file_path_mouse_2 = os.path.join(file_path_1, "example", "collated_pooled_wish_list_test_mouse.txt")

    gene_list, left_primers, right_primers, lib_num = guide.read_gene_list_pooled(file=file_path_5)
    primer_df = pd.DataFrame({'gene_symbol': gene_list,
                              'left_primers': left_primers,
                              'right_primers': right_primers,
                              'lib_num': lib_num})

    glm, lpm, rpm, lnm = guide.read_gene_list_pooled(file=file_path_mouse_2)
    primer_df_mouse = pd.DataFrame({'gene_symbol': glm,
                                    'left_primers': lpm,
                                    'right_primers': rpm,
                                    'lib_num': lnm})

    file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path_2 = os.path.join(file_path_1, "example", "gene_list.txt")
    gene_list = guide.read_gene_list(file_path_2)
    assert gene_list[0] == "STAT3"
    guide.order_guides(gene_list,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="human",
                       primer_df=primer_df)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="mouse",
                       primer_df=primer_df_mouse)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=False,
                       organism="mouse",
                       primer_df=primer_df_mouse)
    guide.order_guides(glm,
                       name="Test",
                       ai_status="a",
                       guides_per_gene=5,
                       order_format="pooled",
                       base_dir=os.path.join(file_path_1, "example"),
                       kampmann_lab=True,
                       organism="mouse",
                       primer_df=primer_df_mouse)
