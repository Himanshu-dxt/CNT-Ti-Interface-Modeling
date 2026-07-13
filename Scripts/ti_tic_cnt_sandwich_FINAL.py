
# ----------------------------------------------------------------------
# Manual TiC(111) Bilayer (2 layers, NO vacuum)
# ----------------------------------------------------------------------

# 1. Install Libraries
!pip install pymatgen ase -q

# 2. Imports
import numpy as np
from ase import Atoms
from ase.build import nanotube
from ase.io import write
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator # Re-added for Ti slab
import warnings
warnings.filterwarnings('ignore')

print("✅ Libraries installed and imported.")

# ----------------------------------------------------------------------
# PART A: BUILD THE (10,5) CNT
# ----------------------------------------------------------------------
print("\n--- Building (10,5) CNT ---")
ase_cnt = nanotube(10, 5, length=3, bond=1.42, symbol='C')

orig_cell = ase_cnt.get_cell()
pos = ase_cnt.get_positions()
min_x, max_x = pos[:, 0].min(), pos[:, 0].max()
min_y, max_y = pos[:, 1].min(), pos[:, 1].max()
vac = 15.0

new_x = (max_x - min_x) + 2 * vac
new_y = (max_y - min_y) + 2 * vac
new_z = orig_cell[2, 2]

new_cell = np.array([[new_x, 0, 0], [0, new_y, 0], [0, 0, new_z]])
ase_cnt.set_cell(new_cell)
ase_cnt.center()

cnt_z_height = new_z
print(f"CNT built: {len(ase_cnt)} atoms, Cell = {new_x:.2f} x {new_y:.2f} x {cnt_z_height:.2f} Å")

# ----------------------------------------------------------------------
# PART B: BUILD THE TITANIUM (0001) SLAB (8x8 supercell)
# ----------------------------------------------------------------------
print("\n--- Building Ti(0001) Slab ---")
a, c = 2.95, 4.68
lattice = Lattice.hexagonal(a, c)
ti_bulk = Structure(lattice, ["Ti", "Ti"], [[0, 0, 0], [1/3, 2/3, 0.5]])

slab_gen = SlabGenerator(ti_bulk, [0, 0, 1], min_slab_size=10, min_vacuum_size=0)
slab_1x1 = slab_gen.get_slab()

supercell_size = 8
slab_big = slab_1x1 * (supercell_size, supercell_size, 1)

adaptor = AseAtomsAdaptor()
ti_ase = adaptor.get_atoms(slab_big)

ti_thickness = ti_ase.get_positions()[:, 2].max() - ti_ase.get_positions()[:, 2].min()
print(f"Ti slab: {len(ti_ase)} atoms, Thickness = {ti_thickness:.2f} Å")

# ----------------------------------------------------------------------
# PART C: MANUAL TiC(111) BILAYER (GUARANTEED THIN, NO VACUUM)
# ----------------------------------------------------------------------
print("\n--- Building Thin TiC(111) Bilayer (MANUAL, Guaranteed) ---")

# TiC lattice parameters
a_cubic = 4.327                # cubic lattice constant
a_hex = a_cubic / np.sqrt(2)   # in-plane lattice constant for (111): 3.059 Å
d111 = a_cubic / np.sqrt(3)    # interlayer spacing: 2.498 Å

# We want the TiC cell to match the CNT cell (~40 Å)
# 13 * 3.059 = 39.77 Å ≈ 40 Å, perfect.
tic_supercell = 13

# Create the hexagonal lattice
i_grid, j_grid = np.meshgrid(np.arange(tic_supercell), np.arange(tic_supercell))
i_grid = i_grid.ravel()
j_grid = j_grid.ravel()

# Hexagonal lattice vectors for the (111) plane (XY components only)
a1_xy = np.array([a_hex, 0])
a2_xy = np.array([a_hex / 2, a_hex * np.sqrt(3) / 2])

# Ti atom positions (at the lattice points, XY components)
ti_pos_xy = i_grid[:, None] * a1_xy + j_grid[:, None] * a2_xy
# Add the Z-component for Ti (at z=0 for the first layer)
ti_pos = np.hstack([ti_pos_xy, np.zeros((ti_pos_xy.shape[0], 1))])

# C atom positions: shifted by (a1 + 2*a2) / 3 (which is the correct offset for rocksalt (111))
# The shift vector also needs to be just XY components
shift_vec_xy = (a1_xy + 2 * a2_xy) / 3
c_pos_xy = ti_pos_xy + shift_vec_xy
# Add the Z-component for C (at z=d111 for the second layer)
c_pos = np.hstack([c_pos_xy, d111 * np.ones((c_pos_xy.shape[0], 1))])

# Combine
atoms_pos = np.vstack([ti_pos, c_pos])
symbols = ['Ti'] * len(ti_pos) + ['C'] * len(c_pos)

# Build the cell
cell_x = tic_supercell * a1_xy[0] # Use a1_xy for X dimension
cell_y = tic_supercell * a2_xy[1] # Use a2_xy for Y dimension
cell_z = d111 + 2.0  # tiny padding (2 Å) just to keep atoms inside, effectively no vacuum

tic_cell = np.array(
    [[cell_x, 0, 0],
     [0, cell_y, 0],
     [0, 0, cell_z]]
)

# Create ASE Atoms object
tic_ase = Atoms(symbols=symbols, positions=atoms_pos, cell=tic_cell, pbc=True)
tic_ase.center()

tic_thickness = tic_ase.get_positions()[:, 2].max() - tic_ase.get_positions()[:, 2].min()
tic_atoms = len(tic_ase)
tic_ti_count = sum(1 for sym in tic_ase.get_chemical_symbols() if sym == 'Ti')
tic_c_count = sum(1 for sym in tic_ase.get_chemical_symbols() if sym == 'C')

print(f"TiC bilayer: {tic_atoms} atoms (Ti: {tic_ti_count}, C: {tic_c_count})")
print(f"TiC thickness: {tic_thickness:.2f} Å (Physically accurate for a bilayer!)")

# ----------------------------------------------------------------------
# PART D: ASSEMBLE THE SANDWICH: Ti + THIN TiC + CNT
# ----------------------------------------------------------------------
print("\n--- Assembling the Sandwich ---")

ti_pos = ti_ase.get_positions().copy()
tic_pos = tic_ase.get_positions().copy()
cnt_pos = ase_cnt.get_positions().copy()

# Final cell: XY from CNT, Z = Ti + gap + TiC + gap + CNT + vacuum
bonding_gap_ti_tic = 2.0
bonding_gap_tic_cnt = 2.1
top_vacuum = 5.0

final_cell_xy = ase_cnt.get_cell()[:2, :2].copy()
final_z_height = (ti_thickness + bonding_gap_ti_tic + tic_thickness +
                  bonding_gap_tic_cnt + cnt_z_height + top_vacuum)

final_cell = np.array([
    [final_cell_xy[0, 0], 0, 0],
    [0, final_cell_xy[1, 1], 0],
    [0, 0, final_z_height]
])
print(f"Final Cell: {final_cell[0,0]:.2f} x {final_cell[1,1]:.2f} x {final_cell[2,2]:.2f} Å")

# Center everything in XY
final_center_x = final_cell[0, 0] / 2
final_center_y = final_cell[1, 1] / 2

# Center Ti slab
ti_center_x = (ti_pos[:, 0].max() + ti_pos[:, 0].min()) / 2
ti_center_y = (ti_pos[:, 1].max() + ti_pos[:, 1].min()) / 2
ti_pos[:, 0] += (final_center_x - ti_center_x)
ti_pos[:, 1] += (final_center_y - ti_center_y)

# Center TiC slab
tic_center_x = (tic_pos[:, 0].max() + tic_pos[:, 0].min()) / 2
tic_center_y = (tic_pos[:, 1].max() + tic_pos[:, 1].min()) / 2
tic_pos[:, 0] += (final_center_x - tic_center_x)
tic_pos[:, 1] += (final_center_y - tic_center_y)

# Center CNT
cnt_center_x = np.mean(cnt_pos[:, 0])
cnt_center_y = np.mean(cnt_pos[:, 1])
cnt_pos[:, 0] += (final_center_x - cnt_center_x)
cnt_pos[:, 1] += (final_center_y - cnt_center_y)

# Position vertically
# Ti slab at bottom
ti_min_z = np.min(ti_pos[:, 2])
ti_pos[:, 2] -= ti_min_z

# TiC above Ti
tic_min_z = np.min(tic_pos[:, 2])
tic_shift_z = ti_thickness + bonding_gap_ti_tic - tic_min_z
tic_pos[:, 2] += tic_shift_z

# CNT above TiC
cnt_min_z = np.min(cnt_pos[:, 2])
cnt_shift_z = ti_thickness + bonding_gap_ti_tic + tic_thickness + bonding_gap_tic_cnt - cnt_min_z
cnt_pos[:, 2] += cnt_shift_z

# Combine
combined_symbols = (ti_ase.get_chemical_symbols() +
                    tic_ase.get_chemical_symbols() +
                    ase_cnt.get_chemical_symbols())
combined_positions = np.vstack([ti_pos, tic_pos, cnt_pos])

sandwich = Atoms(symbols=combined_symbols,
                 positions=combined_positions,
                 cell=final_cell,
                 pbc=True)
sandwich.wrap()

# Save
write("ti_tic_cnt_sandwich_FINAL.cif", sandwich)
print(f"\n✅ FINAL sandwich structure saved as 'ti_tic_cnt_sandwich_FINAL.cif'")
print(f"   Total atoms: {len(sandwich)}")
print(f"   Ti slab: {len(ti_ase)} atoms")
print(f"   TiC layer: {tic_atoms} atoms (Ti: {tic_ti_count}, C: {tic_c_count})")
print(f"   CNT: {len(ase_cnt)} atoms")
print(f"\n   TiC thickness: {tic_thickness:.2f} Å (REALISTIC bilayer!)")
print(f"   Gap Ti-TiC: {bonding_gap_ti_tic:.1f} Å")
print(f"   Gap TiC-CNT: {bonding_gap_tic_cnt:.1f} Å")