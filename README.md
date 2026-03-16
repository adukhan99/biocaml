# BIOCAML

BIOCAML is a command-line tool and OCaml library designed for structural bioinformatics, molecular dynamics (MD) analysis, and protein structure manipulation. It provides a comprehensive suite of utilities for handling common file formats like PDB, CRD, XYZ, and DCD.

## Key Features

- **Multi-Format Support**: Read and write capabilities for PDB, XYZ, and CRD files.
- **Trajectory Analysis**: Parse and analyze CHARMM/NAMD DCD trajectory files.
- **Structural Alignment**: Perform optimal structure alignment using the Kabsch algorithm (via quaternion method).
- **Principal Component Analysis (PCA)**: Analyze conformational changes across multi-frame trajectories.
- **RMSD Calculation**: Compute Root-Mean-Square Deviation between structures or trajectory frames.
- **Geometric Operations**: Calculate bounding boxes, geometric centers, and diameters of molecular systems.
- **Atom Selection**: Powerful filtering system to extract specific chains, residue types, or elements from PDB files.
- **Solvation Tools**: Estimate solvent requirements for molecular systems based on XYZ coordinates and target density.

## Installation

BIOCAML is built using the Dune build system.

### Prerequisites

- OCaml 5.2+
- [Dune](https://dune.build/)
- [Opam](https://opam.ocaml.org/)

### Build from Source

```bash
git clone https://github.com/adukhan99/biocaml.git
cd biocaml
dune build
```

## Usage

The `biocaml` executable provides several subcommands for common tasks.

### PDB Utilities

Summarize a PDB file:
```bash
biocaml pdb-info protein.pdb
```

Calculate geometric bounds of a structure:
```bash
biocaml pdb-bounds complex.pdb
```

Filter a PDB by chain and residue range:
```bash
biocaml pdb-select input.pdb -o output.pdb --chain A --res-min 1 --res-max 100
```

### Trajectory & Analysis

Calculate RMSD between two PDB structures:
```bash
biocaml pdb-rmsd reference.pdb target.pdb
```

Align frames in an XYZ trajectory to a reference:
```bash
biocaml xyz-fit trajectory.xyz -o aligned.xyz --ref reference.xyz
```

Run PCA on a trajectory and output principal components (scores) to CSV:
```bash
biocaml xyz-pca trajectory.xyz --scores scores.csv --modes modes.xyz --k 3
```

### File Conversions

Convert a trajectory from DCD to XYZ format:
```bash
biocaml dcd-to-xyz input.dcd -o output.xyz
```

Convert PDB to CHARMM CRD format:
```bash
biocaml pdb-to-crd molecule.pdb -o molecule.crd
```

## Library Documentation

BIOCAML can also be used as a library in your own OCaml projects. Add `BIOCAML` to your `dune` file dependencies.

```ocaml
open Bio_pdb
open Bio_analysis

let atoms = read_pdb_atoms "molecule.pdb" in
let coords = atoms_to_coords atoms in
match center coords with
| Some c -> Printf.printf "Center: %f %f %f\n" c.x c.y c.z
| None -> ()
```

## License

BIOCAML is released under the [BSD 2-Clause License](LICENSE).