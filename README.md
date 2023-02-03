# Code repository for "Relaxing Hardware Requirements for Surface Code Circuits using Time-dynamics"

This repository contains the code used to generate circuits, and collect sample statistics,
used in the paper "Relaxing Hardware Requirements for Surface Code Circuits using Time-dynamics".

## reproducing figures and statistics

Sitting next to this README, at the root of the repository,
are several scripts:
`step1_generate_circuits.sh`,
`step2_collect_stats.sh`
`step3_derive_fused_stats.sh`
and `step4_plots.sh`.
Most plots from the paper can be reproduced by setting up a python
environment and then running these scripts in order.
**One major difference is that the collection script has been changed to use
pymatching as the decoder, instead of a proprietary internal decoder
that performs correlated decoding.**
The logical error rates of this decoder as slightly better than pymatching.
Also, the collection script has been configured to take a million shots instead
of a hundred million, so that it finishes faster.
Interested users can increase this number.

Assuming a linux-like system with a working python installation, and a bit of luck,
the following script should in principle reproduce the results from the paper after running
for a few days:

```bash
# Preparation: make virtual environment
python -m venv .venv
source .venv/bin/activate
# Preparation: install dependencies
sudo apt install parallel  # scripts use gnu-parallel to distribute circuit generatoin and plotting work.
pip install -r requirements.txt

./step1_generate_circuits.sh
./step2_collect_stats  # This one might take a week!
# Maybe edit step2 to take fewer shots, or edit step1 to make fewer circuits?

./step3_derive_fused_stats.sh
./step4_plot_stats
```

## directory structure

- `.`: top level of repository, with this README and the `step#` scripts
- `./tools`: scripts used by the top-level `step#` scripts to perform tasks such as generating specific circuits and plots
- `./src`: source root of the code; the directory to include in `PYTHONPATH`
- `./src/midout` the python package containing the circuit generation and debugging code (the name is a Silicon Valley reference, due to many of the circuits being built by conceptually starting from the mid cycle state)
- `./src/midout/planar`: location of code producing surface code circuits on a bounded grid
- `./src/midout/toric`: location of code producing surface code circuits on a periodic grid
- `./src/midout/walking`: location of code producing walking surface code circuits
- `./src/midout/gen`, `./src/midout/util`: utility code for generating and debugging circuits
- `./out`: Scripts are configured to create this directory and write their output to various locations within it.

Unfortunately, there isn't really a consistent "style" to the circuit generation code.
Every circuit is done a bit differently, because we were actively learning and experimenting with different
techniques for making them.
This makes the code a bit of a cacophony of concepts, some of which definitely worked better than others.
