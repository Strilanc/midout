from midout.walking import diamonds
from midout.walking.util import Basis, MarkerType


def test_diamond_state_filtering():
    # make some diamonds
    d0x = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.x, marker_type=MarkerType.expanding
    )
    d0z = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.z, marker_type=MarkerType.expanding
    )
    d1x = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.x, marker_type=MarkerType.contracting
    )
    d1z = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.z, marker_type=MarkerType.contracting
    )
    d2x = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j), basis=Basis.x, marker_type=MarkerType.expanding
    )
    d2z = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j), basis=Basis.z, marker_type=MarkerType.expanding
    )
    d3x = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.x, marker_type=MarkerType.error
    )
    d3z = diamonds.Diamond(
        qubits=(0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j), basis=Basis.z, marker_type=MarkerType.error
    )

    all_diamonds = {d0x, d0z, d1x, d1z, d2x, d2z, d3x, d3z}
    x_diamonds = {d0x, d1x, d2x, d3x}
    z_diamonds = {d0z, d1z, d2z, d3z}
    expanding_diamonds = {d0x, d0z, d2x, d2z}
    contracting_diamonds = {d1x, d1z}
    boundary_diamonds = {d2x, d2z}
    error_diamonds = {d3x, d3z}

    ds = diamonds.DiamondState(diamonds=frozenset(all_diamonds), observables={})

    assert ds.get_diamonds() == all_diamonds
    assert ds.get_diamonds(filter_basis=Basis.x) == x_diamonds
    assert ds.get_diamonds(filter_basis=Basis.z) == z_diamonds
    assert ds.get_diamonds(filter_marker=MarkerType.expanding) == expanding_diamonds
    assert ds.get_diamonds(filter_marker=MarkerType.contracting) == contracting_diamonds
    assert ds.get_diamonds(filter_marker=[MarkerType.expanding, MarkerType.error]) == (
        expanding_diamonds | error_diamonds
    )
    assert ds.get_boundary_diamonds() == boundary_diamonds
