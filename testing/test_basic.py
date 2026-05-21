import pandas as pd

from pyguide import guide


def test_reverse_compliment():
    seq1 = "ATTGCA"
    seq2 = "GCATC"
    rc_seq1 = guide.reverse_compliment(seq1)
    rc_seq2 = guide.reverse_compliment(seq2)
    assert rc_seq1 == "TGCAAT"
    assert rc_seq2 == "GATGC"


def test_filter_database_both_empty():
    df = pd.DataFrame({
        'name': ['guide1', 'guide2'],
        'gene': ['STAT3', 'APOE'],
        'score': [1.0, 2.0],
        'seq': ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
    })
    result = guide.filter_database(df, guide_ids=[], gene_names=[])
    assert len(result) == 0
    assert list(result.columns) == list(df.columns)


def test_make_query_map_string_alias():
    query_df = pd.DataFrame(
        {'symbol': ['APOE'], 'alias': ['APOE2']},
        index=['APOE']
    )
    result = guide.make_query_map(query_df)
    assert 'APOE2' in result
    assert result['APOE2'] == 'APOE'
    assert 'A' not in result  # must NOT be split into individual characters
    assert 'P' not in result


def test_filter_cloned_guides_nonsequential_index():
    # DataFrame with labels >> shape[0] simulates post-sort state
    df = pd.DataFrame({
        'gene':  ['STAT3', 'STAT3', 'APOE', 'APOE'],
        'score': [10.0,     8.0,    9.0,    7.0],
        'seq':   ['AAAA',  'TTTT', 'CCCC', 'GGGG'],
        'name':  ['s1',    's2',   'a1',   'a2'],
    }, index=[100, 101, 102, 103])
    cloned = {'STAT3': 1}
    result = guide.filter_cloned_guides(df, cloned)
    stat3_rows = result[result['gene'] == 'STAT3']
    assert len(stat3_rows) == 1
    # Top STAT3 guide (score=10) was already cloned; only score=8 remains
    assert stat3_rows.iloc[0]['score'] == 8.0


def test_split_dataframe_exact_multiple():
    # 96 rows with chunk_size=96 should produce exactly 1 chunk, not 2
    df = pd.DataFrame({'x': range(96)})
    chunks = guide.split_dataframe(df, chunk_size=96)
    assert len(chunks) == 1
    assert len(chunks[0]) == 96


def test_split_dataframe_partial():
    df = pd.DataFrame({'x': range(100)})
    chunks = guide.split_dataframe(df, chunk_size=96)
    assert len(chunks) == 2
    assert len(chunks[0]) == 96
    assert len(chunks[1]) == 4


def test_split_dataframe_empty():
    df = pd.DataFrame({'x': []})
    chunks = guide.split_dataframe(df, chunk_size=96)
    assert len(chunks) == 0

