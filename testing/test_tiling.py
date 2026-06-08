import shutil

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from pyguide.tiling import (
    parse_coordinates,
    reverse_complement,
    fetch_sequence_ucsc,
    find_ngg_guides,
    run_guidescan,
    apply_filters,
)


def test_import_tiling():
    from pyguide import tiling
    assert hasattr(tiling, 'parse_coordinates')


def test_parse_coordinates_valid():
    chrom, start, end = parse_coordinates("chr1:12345-12395")
    assert chrom == "chr1"
    assert start == 12345
    assert end == 12395


def test_parse_coordinates_invalid_format():
    with pytest.raises(ValueError, match="Expected format"):
        parse_coordinates("12345-12395")


def test_parse_coordinates_invalid_range():
    with pytest.raises(ValueError, match="start must be less than end"):
        parse_coordinates("chr1:12395-12345")


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCTA") == "TAGC"


def test_find_ngg_guides_forward():
    # Sequence with a single NGG at position 20-22 (0-indexed)
    # spacer = first 20 nt, PAM = AGG
    seq = "ACGTACGTACGTACGTACGTAGG"  # 23nt: 20nt spacer + AGG
    guides = find_ngg_guides("chr1", 1, seq)
    fwd = [g for g in guides if g['match_strand'] == '+']
    assert len(fwd) >= 1
    assert fwd[0]['sequence'] == "ACGTACGTACGTACGTACGT"
    assert fwd[0]['match_strand'] == '+'
    assert fwd[0]['match_position'] == 1


def test_find_ngg_guides_reverse():
    # CCN on + strand means NGG guide on - strand
    # CC at positions 0-1 means spacer on - strand = revcomp(seq[3:23])
    seq = "CCT" + "ACGTACGTACGTACGTACGT"  # 23nt: CCN + 20nt
    guides = find_ngg_guides("chr1", 1, seq)
    rev = [g for g in guides if g['match_strand'] == '-']
    assert len(rev) >= 1
    expected_spacer = reverse_complement("ACGTACGTACGTACGTACGT")
    assert rev[0]['sequence'] == expected_spacer
    assert rev[0]['match_strand'] == '-'


def test_find_ngg_guides_excludes_n():
    # Spacer with N in it should be excluded
    seq = "ACGTACGTACGTACGNACGTAGG"
    guides = find_ngg_guides("chr1", 1, seq)
    assert not any('N' in g['sequence'] for g in guides)


def test_run_guidescan_missing_binary():
    guides = [
        {'sequence': 'ACGTACGTACGTACGTACGT', 'match_chrm': 'chr1',
         'match_position': 100, 'match_strand': '+'},
    ]
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="guidescan.*not found"):
            run_guidescan(guides, index_path="/fake/index")


def test_run_guidescan_returns_dataframe():
    """run_guidescan merges specificity from CSV onto guide list."""
    guides = [
        {'sequence': 'ACGTACGTACGTACGTACGT', 'match_chrm': 'chr1',
         'match_position': 100, 'match_strand': '+'},
    ]
    fake_csv = (
        "id,sequence,match_chrm,match_position,match_strand,match_distance,specificity\n"
        "guide_0,ACGTACGTACGTACGTACGT,chr1,99,+,0,0.85\n"
    )

    def fake_subprocess(cmd, **kwargs):
        # Write what guidescan would write to --output
        out_idx = cmd.index('--output')
        with open(cmd[out_idx + 1], 'w') as f:
            f.write(fake_csv)
        m = MagicMock()
        m.returncode = 0
        m.stderr = ''
        return m

    with patch('shutil.which', return_value='/usr/bin/guidescan'):
        with patch('subprocess.run', side_effect=fake_subprocess):
            df = run_guidescan(guides, index_path="/fake/index")
    assert 'specificity' in df.columns
    assert len(df) == 1
    assert df.iloc[0]['specificity'] == pytest.approx(0.85)


def _make_guide_df(sequences, specificities=None):
    if specificities is None:
        specificities = [0.9] * len(sequences)
    return pd.DataFrame({
        'sequence': sequences,
        'match_chrm': ['chr1'] * len(sequences),
        'match_position': list(range(len(sequences))),
        'match_strand': ['+'] * len(sequences),
        'specificity': specificities,
    })


def test_apply_filters_removes_tttt():
    df = _make_guide_df(['ACGTTTTTACGTACGTACGT', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_bstxi():
    df = _make_guide_df(['CCACCTTGTTGACGTACGTA', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_bpi1102i():
    df = _make_guide_df(['GTTTAAGAGCTAAGCTGGAC', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_low_specificity():
    df = _make_guide_df(
        ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
        specificities=[0.1, 0.9],
    )
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'TGCATGCATGCATGCATGCA'


def test_apply_filters_removes_nan_specificity():
    df = _make_guide_df(['ACGTACGTACGTACGTACGT'], specificities=[float('nan')])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 0
