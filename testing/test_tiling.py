import pytest
import pandas as pd
from pyguide.tiling import (
    parse_coordinates,
    reverse_complement,
    fetch_sequence_ucsc,
    find_ngg_guides,
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
