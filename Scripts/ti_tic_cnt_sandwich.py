# ----------------------------------------------------------------------
# CNT(10,5) on TiC(111) on Ti(0001) Interface
# Builds the full experimentally-relevant sandwich structure
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
from pymatgen.core.surface import SlabGenerator
import warnings
warnings.filterwarnings('ignore')

print("✅ Libraries installed and imported.")

# ----------------------------------------------------------------------
# PART A: BUILD THE (10,5) CNT (with vacuum padding) - SAME AS BEFORE
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
# PART B: BUILD THE TITANIUM (0001) SLAB (8x8 supercell) - SAME AS BEFORE
# ----------------------------------------------------------------------
print("\n--- Building Ti(0001) Slab ---")
a, c = 2.95, 4.68
lattice = Lattice.hexagonal(a, c)
ti_bulk = Structure(lattice, ["Ti", "Ti"], [[0, 0, 0], [1/3, 2/3, 0.5]])

slab_gen = SlabGenerator(ti_bulk, [0, 0, 1], min_slab_size=10, min_vacuum_size=0)
slab_1x1 = slab_gen.get_slab()

# Expand to 8x8 supercell to match CNT size
supercell_size = 8
slab_big = slab_1x1 * (supercell_size, supercell_size, 1)

adaptor = AseAtomsAdaptor()
ti_ase = adaptor.get_atoms(slab_big)

ti_thickness = ti_ase.get_positions()[:, 2].max() - ti_ase.get_positions()[:, 2].min()
print(f"Ti slab: {len(ti_ase)} atoms, Thickness = {ti_thickness:.2f} Å")

# ----------------------------------------------------------------------
# PART C: BUILD THE TiC(111) LAYER  (NEW!)
# ----------------------------------------------------------------------
print("\n--- Building TiC(111) Layer ---")

# TiC has rock-salt (NaCl) structure, space group Fm-3m
# Lattice parameter: 4.327 Å at ambient conditions [7†L28-L30]
tic_lattice_param = 4.327

# Create TiC bulk structure (rock-salt)
# Ti at (0,0,0) and C at (0.5, 0.5, 0.5)
tic_bulk = Structure(
    Lattice.cubic(tic_lattice_param),
    ["Ti", "C"],
    [[0, 0, 0], [0.5, 0.5, 0.5]]
)

# Generate (111) surface slab
# Note: TiC(111) is polar - it has two different surface terminations (Ti-terminated and C-terminated) [1†L4-L7]
# We'll generate a slab with both terminations
tic_slab_gen = SlabGenerator(
    tic_bulk,
    [1, 1, 1],          # Miller index for (111) surface
    min_slab_size=8,    # Minimum slab thickness in Angstroms
    min_vacuum_size=0   # No vacuum yet - we'll add it globally
)

# Get the most stable slab (usually the one with both surfaces terminated the same way)
tic_slabs = tic_slab_gen.get_slabs()
tic_slab = tic_slabs[0]  # Take the first (most symmetric) slab

# Expand TiC slab to match the CNT's XY dimensions
# TiC lattice constant is ~4.327 Å, CNT cell is ~40 Å
# So we need roughly 40/4.327 ≈ 9 unit cells
tic_supercell_size = 9
tic_slab_big = tic_slab * (tic_supercell_size, tic_supercell_size, 1)

tic_ase = adaptor.get_atoms(tic_slab_big)
tic_thickness = tic_ase.get_positions()[:, 2].max() - tic_ase.get_positions()[:, 2].min()
print(f"TiC slab: {len(tic_ase)} atoms, Thickness = {tic_thickness:.2f} Å")

# ----------------------------------------------------------------------
# PART D: ASSEMBLE THE SANDWICH: Ti + TiC + CNT
# ----------------------------------------------------------------------
print("\n--- Assembling the Sandwich: Ti + TiC + CNT ---")

# 1. Get positions from all three components
ti_pos = ti_ase.get_positions().copy()
tic_pos = tic_ase.get_positions().copy()
cnt_pos = ase_cnt.get_positions().copy()

# 2. Define the final cell dimensions
# XY: Use the CNT's XY cell (which has vacuum padding)
final_cell_xy = ase_cnt.get_cell()[:2, :2].copy()

# Z: Ti slab + TiC slab + bonding gap (Ti-TiC) + CNT + top vacuum
# Bonding gap between Ti and TiC: ~2.0-2.5 Å (typical metal-ceramic interface)
bonding_gap_ti_tic = 2.0
# Bonding gap between TiC and CNT: ~2.1 Å (same as before)
bonding_gap_tic_cnt = 2.1
# Top vacuum above CNT
top_vacuum = 5.0

final_z_height = (ti_thickness + bonding_gap_ti_tic + tic_thickness + 
                  bonding_gap_tic_cnt + cnt_z_height + top_vacuum)

final_cell = np.array([
    [final_cell_xy[0, 0], 0, 0],
    [0, final_cell_xy[1, 1], 0],
    [0, 0, final_z_height]
])
print(f"Final Cell: {final_cell[0,0]:.2f} x {final_cell[1,1]:.2f} x {final_cell[2,2]:.2f} Å")

# 3. Center everything in XY
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

# 4. Position vertically (Z-axis)
# Ti slab goes at the bottom (Z=0)
# Shift Ti slab so its bottom is at Z=0
ti_min_z = np.min(ti_pos[:, 2])
ti_pos[:, 2] -= ti_min_z

# TiC slab goes above Ti with a gap
tic_min_z = np.min(tic_pos[:, 2])
tic_shift_z = ti_thickness + bonding_gap_ti_tic - tic_min_z
tic_pos[:, 2] += tic_shift_z

# CNT goes above TiC with a gap
cnt_min_z = np.min(cnt_pos[:, 2])
cnt_shift_z = ti_thickness + bonding_gap_ti_tic + tic_thickness + bonding_gap_tic_cnt - cnt_min_z
cnt_pos[:, 2] += cnt_shift_z

# 5. Combine everything
combined_symbols = (ti_ase.get_chemical_symbols() + 
                    tic_ase.get_chemical_symbols() + 
                    ase_cnt.get_chemical_symbols())
combined_positions = np.vstack([ti_pos, tic_pos, cnt_pos])

sandwich = Atoms(symbols=combined_symbols,
                 positions=combined_positions,
                 cell=final_cell,
                 pbc=True)
sandwich.wrap()

# Save the structure
write("ti_tic_cnt_sandwich.cif", sandwich)
print(f"\n✅ Sandwich structure saved as 'ti_tic_cnt_sandwich.cif'")
print(f"   Total atoms: {len(sandwich)}")
print(f"   Ti atoms: {len(ti_ase)}")
print(f"   TiC atoms: {len(tic_ase)} (Ti + C)")
print(f"   C atoms (CNT): {len(ase_cnt)}")
print(f"\n   Structure: Ti({len(ti_ase)}) + TiC({len(tic_ase)}) + CNT({len(ase_cnt)})")
print(f"   Gap Ti-TiC: {bonding_gap_ti_tic:.1f} Å")
print(f"   Gap TiC-CNT: {bonding_gap_tic_cnt:.1f} Å")