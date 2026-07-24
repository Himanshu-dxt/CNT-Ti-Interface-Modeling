# Atomistic Modeling of Ti-CNT Interfaces for CNTFET Contacts

**Author:** Himanshu Kumar Dixit  
**Affiliation:** Government Engineering College, Ajmer, India  
**Date:** 13th July, 2026

---

## 📖 Overview
This repository contains a suite of atomistic models for studying metal contacts to semiconducting carbon nanotubes (CNTs). The work addresses the critical challenge of contact resistance in Carbon Nanotube Field-Effect Transistors (CNTFETs) by modeling the interface between a (10,5) chiral CNT and Titanium (Ti) based electrodes. It includes both a simple Ti-CNT interface and a physically realistic Ti-TiC-CNT sandwich structure that accounts for the formation of Titanium Carbide (TiC) upon thermal annealing—a key experimental step for achieving low-resistance ohmic contacts.

The primary objectives are to provide validated, publication-ready geometries for Density Functional Theory (DFT) validation and to demonstrate a reproducible, open-source computational workflow for interface construction.

---

## 🧬 Models Included

| Model | Description | Atoms | File |
| **Simple Ti-CNT Interface** | CNT(10,5) on Ti(0001) surface at 2.1 Å gap. | 548 | `ti_cnt_10_5_interface.cif` |
| **TiC Mediated Interface (Thick)** | CNT on TiC(111) (thick slab) on Ti. | 1196 | `ti_tic_cnt_sandwich.cif` |
| **TiC Mediated Interface (Thin)** | CNT on TiC(111) (thin slab) on Ti. | 710 | `ti_tic_cnt_sandwich_THIN.cif` |
| **TiC Mediated Interface (FINAL)** | CNT on TiC(111) **bilayer** (2.50 Å) on Ti. | **886** | `ti_tic_cnt_sandwich_FINAL.cif` |

### 🔬 The Significance of the FINAL Model
- It captures the **experimentally relevant structure** where TiC forms as a thin bilayer during device fabrication.
- The 2.50 Å thickness matches the theoretical interlayer spacing of TiC(111).
- This structure explains the low contact resistance (< 1.5 kΩ·µm) observed in annealed Ti contacts, despite the work function mismatch (~4.3 eV Ti vs ~4.8 eV CNT).

---

## 📁 Repository Structure
```
CNT-Ti-Interface-Modeling/
│
├── README.md                           # This file
├── requirements.txt                    # Python dependencies
│
├── scripts/                            # Segmented Python scripts for reproducibility
│   ├── 1_ti_cnt_10_5_interface.py      # Builds the simple Ti-CNT interface
│   ├── 2_ti_tic_cnt_sandwich.py        # Builds the thick TiC model
│   ├── 3_ti_tic_cnt_sandwich_THIN.py   # Builds the thin TiC model
│   └── 4_ti_tic_cnt_sandwich_FINAL.py  # Builds the final bilayer sandwich (886 atoms)
│
├── structures/                         # Final .cif files for visualization
│   ├── ti_cnt_10_5_interface.cif
│   ├── ti_tic_cnt_sandwich.cif
│   ├── ti_tic_cnt_sandwich_THIN.cif
│   └── ti_tic_cnt_sandwich_FINAL.cif
│
└── docs/
    └── images/                         # VESTA screenshots of the models
        ├── 10,5_CNT.png
        ├── CNT_Ti_interface_front_view.png
        ├── CNT_Ti_interface_side_view.png
        ├── Ti_slab_front_view.png
        ├── Ti_slab_side_view.png
        ├── Ti_slab_top_view.png
        ├── ti_tic_cnt_sandwich_FINAL_side_view.png
        ├── ti_tic_cnt_sandwich_THIN_side_view.png
        └── ti_tic_cnt_sandwich_side_view.png
```

---

## 🛠️ Methodology

The models are constructed using open-source tools:
- **ASE (Atomic Simulation Environment):** For generating the (10,5) nanotube and handling atomic structures.
- **Pymatgen:** For bulk crystal generation and slab construction.
- **Manual Crystallography:** The TiC(111) bilayer is built atom-by-atom to ensure precise orientation and stoichiometry (1:1 Ti:C).

### Key Parameters
| Parameter | Value |
| :--- | :--- |
| CNT Chirality | (10,5) Semiconducting |
| CNT Diameter | ~1.04 nm |
| CNT Length | 3 Unit Cells (33.81 Å) |
| Ti Slab | 8×8 Ti(0001) Supercell (128 atoms) |
| TiC Orientation | (111) Bilayer (Ti-rich termination) |
| Ti-TiC Gap | 2.0 Å |
| TiC-CNT Gap | 2.1 Å |
| TiC Thickness (Final) | 2.50 Å |
| Final Cell Size | 40.36 × 40.35 × 50.09 Å³ |

---

## 🚀 Getting Started

### Prerequisites
- Python 3.8+
- Required libraries: `ase`, `pymatgen`, `numpy`

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/Himanshu-dxt/CNT-Ti-Interface-Modeling.git
   cd CNT-Ti-Interface-Modeling
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run any of the scripts in the `scripts/` folder, e.g.:
   ```bash
   python scripts/4_ti_tic_cnt_sandwich_FINAL.py
   ```

### Visualization
The resulting `.cif` files can be opened in:
- **VESTA** (Recommended)
- **OVITO**
- **PyMOL**

---

## 📊 Results & Validation
- **Structural Integrity:** Visual inspection confirms the CNT is intact, the TiC bilayer is continuous, and both interfacial gaps are correctly maintained.
- **Physical Consistency:** The TiC bilayer thickness (2.50 Å) matches the theoretical interlayer spacing for TiC(111).
- **Reproducibility:** All structures are generated via self-contained Python scripts without external dependencies on pre-existing data files.

The final model provides a validated foundation for DFT studies on Schottky barrier heights, charge transfer, and the role of interfacial chemistry in contact resistance.

---

## 🔬 Research Context

### The Problem
Titanium has a work function (~4.3 eV) that is significantly lower than the valence band edge of semiconducting CNTs (~4.8 eV). A simple Schottky-Mott model predicts a large barrier to hole injection, suggesting Ti should be a poor contact metal. However, experiments show the opposite: annealed Ti contacts exhibit exceptionally low resistance.

### The Solution
Thermal annealing induces the formation of a thin TiC layer at the interface. TiC is metallic and electronically bridges the Ti and CNT, eliminating the Schottky barrier. This repository models that exact structure, allowing computational validation of the experimental observations.

---

## 📝 Author
**Himanshu Kumar Dixit**  
Mechanical Engineering  
Government Engineering College, Ajmer, India  
📧 24me19@ecajmer.ac.in

---

## 📄 License
This project is open-source and available for academic and research purposes under the MIT License.

---

## 🙏 Acknowledgments
This work was conducted independently as a student project at Government Engineering College, Ajmer. The author gratefully acknowledges the open-source community for developing and maintaining the ASE and Pymatgen libraries.
