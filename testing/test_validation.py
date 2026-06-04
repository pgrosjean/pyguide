import pytest
import os
from pyguide import guide


def test_order_guides_too_many_per_gene():
    with pytest.raises(AssertionError):
        guide.order_guides(
            [], ['STAT3'], name='Test', ai_status='i',
            guides_per_gene=11, order_format='single',
            base_dir='/tmp', check_db=False, organism='human',
        )


def test_order_guides_both_inputs_empty():
    with pytest.raises(AssertionError):
        guide.order_guides(
            [], [], name='Test', ai_status='i',
            guides_per_gene=5, order_format='single',
            base_dir='/tmp', check_db=False, organism='human',
        )


def test_order_guides_invalid_format():
    file_path = os.path.dirname(os.path.abspath(__file__))
    gene_list = ['STAT3']
    with pytest.raises(AssertionError):
        guide.order_guides(
            [], gene_list, name='Test', ai_status='i',
            guides_per_gene=5, order_format='invalid_format',
            base_dir=os.path.join(file_path, 'example'),
            check_db=False, organism='human',
        )
