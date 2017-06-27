
import pytest

from dtocean_electrical.input_utils.utils import set_burial_from_bpi


@pytest.mark.parametrize("test_input, expected", [
    ('hard glacial till', 0),
    ('cemented', 0),
    ('soft rock coral', 0),
    ('hard rock', 0),
    ('gravel cobble', 0),
    ('loose sand', 0.5),
    ('medium sand', 0.5),
    ('dense sand', 0.5),
    ('very soft clay', 1.0),
    ('soft clay', 1.0),
    ('firm clay', 1.0),
    ('stiff clay', 1.0)
    ])
def test_set_burial_from_bpi(test_input, expected):
    
    row = {'layer 1 type': test_input}
    bpi = set_burial_from_bpi(row)
    
    assert bpi == expected


def test_set_burial_from_bpi_error():
    
    row = {'layer 1 type': "shingle"}
    
    with pytest.raises(ValueError):
        set_burial_from_bpi(row)
