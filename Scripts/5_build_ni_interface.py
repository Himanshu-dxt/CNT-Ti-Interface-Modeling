# ----------------------------------------------------------------------
# BUILD Ni-NiC-CNT INTERFACE (Structural Analog to Ti-TiC-CNT)
# Uses ASE and Pymatgen. Runs completely in Google Colab.
# No distortion: CNT built in its own cell, then placed on matching Ni slab.
# Author: Himanshu Dixit
# Date: July 2026
# ----------------------------------------------------------------------

# 1. Install required libraries
!pip install ase pymatgen -q

# 2. Imports
import numpy as np
from ase import Atoms
from ase.build import nanotube, fcc111
from ase.io import write
import warnings
warnings.filterwarnings('ignore')

print("✅ Libraries installed and imported.")

# ----------------------------------------------------------------------
# STEP 1: Build the (10,5) CNT with its natural cell
# ----------------------------------------------------------------------
print("\n--- Building (10,5) CNT ---")
cnt = nanotube(10, 5, length=3, bond=1.42, symbol='C')
cnt_cell = cnt.get_cell().copy()
cnt_pos = cnt.get_positions().copy()
print(f"CNT: {len(cnt)} atoms, cell: {cnt_cell[0,0]:.2f} x {cnt_cell[1,1]:.2f} Å")

# ----------------------------------------------------------------------
# STEP 2: Build the Ni(111) slab (FCC)
# ----------------------------------------------------------------------
print("\n--- Building Ni(111) slab ---")
supercell_size = 12
ni_slab = fcc111('Ni', size=(supercell_size, supercell_size, 2), a=3.52, vacuum=0.0)
ni_pos = ni_slab.get_positions().copy()
ni_top = ni_pos[:, 2].max()
print(f"Ni slab: {len(ni_slab)} atoms, Ni top: {ni_top:.2f} Å")

# ----------------------------------------------------------------------
# STEP 3: Build NiC(111) bilayer (rock-salt)
# ----------------------------------------------------------------------
print("\n--- Building NiC(111) bilayer ---")
a_nic = 4.2 / np.sqrt(2)   # in-plane: 2.97 Å
d_nic = 4.2 / np.sqrt(3)   # interlayer spacing: 2.42 Å
n_nic = 14

nic_positions = []
nic_symbols = []
for i in range(n_nic):
    for j in range(n_nic):
        x = i * a_nic + j * a_nic/2
        y = j * a_nic * np.sqrt(3)/2
        # Ni layer (bottom of NiC)
        nic_positions.append([x, y, 0.0])
        nic_symbols.append('Ni')
        # C layer (top of NiC)
        nic_positions.append([x + a_nic/3, y + a_nic*np.sqrt(3)/6, d_nic])
        nic_symbols.append('C')

nic_pos = np.array(nic_positions)

# Center NiC on the Ni slab
ni_center_xy = np.array([ni_slab.cell[0,0]/2, ni_slab.cell[1,1]/2])
nic_center_xy = np.mean(nic_pos[:, :2], axis=0)
nic_pos[:, :2] += (ni_center_xy - nic_center_xy)

# Place NiC above the Ni slab (2.0 Å gap)
nic_pos[:, 2] += ni_top + 2.0

# ----------------------------------------------------------------------
# STEP 4: Calculate NiC top and place CNT above it
# ----------------------------------------------------------------------
nic_top = nic_pos[:, 2].max()
print(f"NiC top: {nic_top:.2f} Å")

# Center CNT horizontally
cnt_center_xy = np.mean(cnt_pos[:, :2], axis=0)
cnt_pos[:, :2] += (ni_center_xy - cnt_center_xy)

# Place CNT ABOVE the NiC layer (2.1 Å gap)
cnt_bottom = cnt_pos[:, 2].min()
cnt_pos[:, 2] += (nic_top + 2.1 - cnt_bottom)

# ----------------------------------------------------------------------
# STEP 5: Assemble the full interface
# ----------------------------------------------------------------------
print("\n--- Assembling final interface ---")
all_pos = np.vstack([ni_pos, nic_pos, cnt_pos])
all_symbols = ['Ni'] * len(ni_pos) + nic_symbols + ['C'] * len(cnt)

interface = Atoms(symbols=all_symbols, positions=all_pos, cell=ni_slab.cell, pbc=True)
interface.center(vacuum=15.0, axis=2)

# ----------------------------------------------------------------------
# STEP 6: Save and download
# ----------------------------------------------------------------------
write('interface_Ni_structure.cif', interface)
print(f"\n✅ Ni-NiC-CNT interface built (CORRECTED)!")
print(f"   Total atoms: {len(interface)}")
print(f"   Cell dimensions: {interface.cell[0,0]:.2f} x {interface.cell[1,1]:.2f} x {interface.cell[2,2]:.2f} Å")
print("   Saved to interface_Ni_structural_corrected.cif")

print("✅ Structure saved as interface_Ni_structure.cif")
print("\n📥 File downloaded to your local machine.")
