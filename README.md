# CNT(10,5) on Ti(0001) Interface: An Atomistic Model

**Author:** Himanshu Kumar Dixit
**Date:** 13th July, 2026

## Overview
This repository contains the code and structures for building an atomistic model of a (10,5) semiconducting carbon nanotube (CNT) interfaced with a Ti(0001) surface. This is a proof-of-concept for a computational study on metal-CNT contacts for next-generation field-effect transistors (CNTFETs).

The workflow leverages open-source tools: `pymatgen` for structure manipulation[reference:2], `ase` for nanotube generation, and the `Materials Project` API for reference data[reference:3].

## Key Features
- **Chiral CNT Generation:** Scripts to build a (10,5) nanotube with a ~1.04 nm diameter.
- **Realistic Metal Contact:** Construction of a large 8x8 Ti(0001) supercell (128 atoms).
- **Physically Accurate Interface:** The CNT is placed at the correct bonding distance (2.1 Å) above the Ti slab.
- **Open and Reproducible:** All code is provided in Python notebooks and scripts.

## File Structure
- `structures/`: Final `.cif` files of the CNT, Ti slab, and the assembled interface.
- `docs/images/`: Visualization of the final structure (e.g., VESTA screenshots).

## Getting Started
1.  Clone this repository.
2.  Install the required Python packages: `pip install -r requirements.txt`
3.  Run the `notebooks/cnt_ti_interface_builder.ipynb` notebook.

## Results
The final model consists of a **548-atom** system (128 Ti + 420 C). The structure is fully periodic and ready for further computational analysis, such as DFT or ML-based geometry optimization. See `docs/images/` for visualizations.

## Next Steps
- Perform geometry optimization using a suitable force field.
- Validate binding energy with Density Functional Theory (DFT).
- Extend the model to study other contact metals (e.g., Pd, Au).

## Acknowledgments
This work was conducted independently as a student project at Government Engineering College, Ajmer.
