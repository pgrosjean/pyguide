from pyguide import guide

def test_reverse_compliment():
    seq1 = "ATTGCA"
    seq2 = "GCATC"
    rc_seq1 = guide.reverse_compliment(seq1)
    rc_seq2 = guide.reverse_compliment(seq2)
    assert rc_seq1 == "TGCAAT"
    assert rc_seq2 == "GATGC"

