# Homodesmotic Family Reactions
**HFR** uses **Integer Linear Programming** along with **RDKit** and **PuLP** to automate the construction of homodesmotic family reactions. The homodesmotic family reactions can be coupled with **AaronTools** for Gaussian computations. The program is written in Python and designed to run in the command-line terminal. 

The core functionality is implemented in `core.py`. Both **HFR** and **MHFR** are interfaces to `core.py`, providing different input and output options:

- **HFR** offers `view`, `write`, and `count` commands to:
  - Construct reactions
  - Write Gaussian `.com` input files into a designated outfolder
  - Count variables for a specified homodesmotic family reaction level

- **MHFR** processes an input file containing integer levels (1-4) and SMILES strings, then writes the corresponding reactions to an output file.

**Additional Tools**

  - **compute_fe** works within the **HFR** `write` outfolder directory and extracts free energy calculations to determine the **HFR** `input_mol` free energy.

  - **read.py** extracts and prints energies of input Gaussian `.log` file.
