from midout import gen


def test_mul():
    a = 'IIIIXXXXYYYYZZZZ'
    b = 'IXYZ' * 4
    c = 'IXYZXIZYYZIXZYXI'
    a = gen.PauliString({q: p for q, p in enumerate(a) if p != 'I'})
    b = gen.PauliString({q: p for q, p in enumerate(b) if p != 'I'})
    c = gen.PauliString({q: p for q, p in enumerate(c) if p != 'I'})
    assert a * b == c

