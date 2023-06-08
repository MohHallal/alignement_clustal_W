# Sequence Alignment and Phylogenetic Tree Construction

This repository contains Python scripts for performing sequence alignment and constructing phylogenetic trees. The scripts utilize various algorithms such as Needleman-Wunsch, UPGMA, and Neighbor Joining.

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Scripts](#scripts)
- [License](#license)

## Overview

The scripts in this repository provide functionality for sequence alignment and phylogenetic tree construction using different algorithms. The main steps involved are:

1. Reading sequences from a FASTA file.
2. Computing sequence alignment scores using the Needleman-Wunsch algorithm.
3. Building a guide tree using the UPGMA algorithm.
4. Performing global sequence alignment based on the guide tree.
5. Identifying conserved positions in the aligned sequences.
6. Constructing a distance matrix based on the conserved positions.
7. Generating a phylogenetic tree using the Neighbor Joining algorithm.

## Dependencies

The following dependencies are required to run the scripts:

- Python 3.x
- [numpy](https://numpy.org/)
- [biopython](https://biopython.org/)

Please make sure you have these dependencies installed before running the scripts.

## Usage

To use the scripts, follow these steps:

1. Ensure that the required dependencies are installed.
2. Place your sequences in a FASTA file named `opsines` in the same directory as the scripts.
3. Run the `main.py` script using Python: `python main.py`.
4. The script will display the results in the console, including the scores matrix, the guide tree, the multiple alignment, the conserved positions, the distance matrix, and the final phylogenetic tree.

## Scripts

The repository includes the following scripts:

- `alignement_global.py`: Contains functions related to global sequence alignment.
- `alignement_local.py`: Contains functions related to local sequence alignment.
- `neibour_joining.py`: Contains functions for the Neighbor Joining algorithm.
- `upgma.py`: Contains functions for the UPGMA algorithm.
- `main.py`: The main script that orchestrates the sequence alignment and phylogenetic tree construction process.

## License

This project is licensed under the [MIT License](LICENSE).

