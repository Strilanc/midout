import itertools
import os

from midout.walking.circuit import Circuit
from midout.walking.util import Basis, DOWN_RIGHT, UP_LEFT, UP_RIGHT, DOWN_LEFT

CIRCUITS_PATH = "./circuits/"


def make_circuit_files(overwrite=False):
    if not os.path.exists(CIRCUITS_PATH):
        print("MAKING CIRCUIT FOLDER")
        os.makedirs(CIRCUITS_PATH)

    distances = [5, 9, 9, 11, 13, 15, 17, 19, 21]
    error_probs = [None, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.013, 0.015, 0.02]
    types = {
        'standard': [None, None],
        'gliding': [DOWN_RIGHT, DOWN_RIGHT],
        'sliding': [DOWN_RIGHT, UP_RIGHT],
        'wiggling': [DOWN_RIGHT, UP_LEFT],
    }

    for b, d, p, t in itertools.product(Basis, distances, error_probs, types):
        c_filename = f'd={d},p={p},b={b.value},type={t}'
        if not os.path.exists(CIRCUITS_PATH + c_filename) or overwrite:
            print('MAKING CIRCUIT: ' + c_filename)
            
            rounds = types[t] * (d//2)
            if d%2:
                rounds += [None]

            c = Circuit(distance=d, rounds=rounds, basis=b, duids=False).build_stim_circuit(
                errors=p
            )
            if p is not None:
                err_len = len(c.shortest_graphlike_error())
                if err_len != d:
                    raise ValueError(f"PROBLEM: d={d},p={p},b={b.value},type={t} err_len={err_len}")

            with open(CIRCUITS_PATH + c_filename+'.stim', 'w') as f:
                print(c, file=f)

            with open(CIRCUITS_PATH + c_filename+'.svg', 'w') as f:
                print(c.diagram(type="detector-slice-svg", tick=range(1, d*6+2)), file=f)


if __name__ == '__main__':
    make_circuit_files(overwrite=True)
