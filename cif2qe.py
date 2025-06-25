import streamlit as st
from ase.io import read
from ase.data import atomic_masses, chemical_symbols
import tempfile
import os

def ase_symbol_to_z(symbol):
    return chemical_symbols.index(symbol)

def generate_qe_input(atoms, pseudo_map, calc_type, ecutwfc, ecutrho, kpoints):
    species = sorted(set(atoms.get_chemical_symbols()))
    missing = [el for el in species if el not in pseudo_map]
    if missing:
        st.error(f"Missing pseudopotentials for: {', '.join(missing)}")
        return None, None

    lines = []
    lines.append("&control")
    lines.append(f"  calculation = '{calc_type}',")
    lines.append("  prefix = 'calc',")
    lines.append("  outdir = './tmp',")
    lines.append("/\n")

    lines.append("&system")
    lines.append(f"  ibrav = 0,")
    lines.append(f"  nat = {len(atoms)},")
    lines.append(f"  ntyp = {len(species)},")
    lines.append(f"  ecutwfc = {ecutwfc},")
    lines.append(f"  ecutrho = {ecutrho},")
    lines.append("/\n")

    lines.append("&electrons")
    lines.append("  conv_thr = 1.0d-6")
    lines.append("/\n")

    lines.append("CELL_PARAMETERS angstrom")
    for vec in atoms.get_cell():
        lines.append("  {:.10f} {:.10f} {:.10f}".format(*vec))

    lines.append("\nATOMIC_SPECIES")
    for el in species:
        z = ase_symbol_to_z(el)
        mass = atomic_masses[z]
        lines.append(f"{el}  {mass:.4f}  {pseudo_map[el].name}")

    lines.append("\nATOMIC_POSITIONS angstrom")
    for atom in atoms:
        pos = atom.position
        lines.append(f"{atom.symbol}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")

    lines.append("\nK_POINTS automatic")
    lines.append("  " + " ".join(kpoints) + " 0 0 0")

    qe_input_str = "\n".join(lines)
    run_script = f"""#!/bin/bash
# Run Quantum ESPRESSO
pw.x < qe_input.in > output.log
"""

    return qe_input_str, run_script

# ---------- Streamlit UI ----------
st.title("🧪 CIF to Quantum ESPRESSO Input Generator")
st.markdown("Upload a CIF file, assign pseudopotentials for each atom, and generate the QE input file.")

# Upload CIF file
cif_file = st.file_uploader("📁 Upload CIF file", type=["cif"])

if cif_file:
    try:
        atoms = read(cif_file)
        species = sorted(set(atoms.get_chemical_symbols()))
        st.success(f"✅ Parsed CIF. Detected elements: {', '.join(species)}")
    except Exception as e:
        st.error(f"❌ Failed to read CIF file: {str(e)}")
        atoms = None
else:
    atoms = None
    species = []

# Upload UPF for each atom type
pseudo_map = {}
if atoms:
    st.subheader("🔗 Upload Pseudopotentials")
    for el in species:
        pseudo = st.file_uploader(f"Select pseudopotential for `{el}`", type=["UPF", "upf"], key=el)
        if pseudo:
            pseudo_map[el] = pseudo

# QE Parameters
st.subheader("⚙ QE Input Parameters")
col1, col2, col3 = st.columns(3)
with col1:
    calc_type = st.selectbox("Calculation Type", ["scf", "relax", "vc-relax", "nscf"])
with col2:
    ecutwfc = st.number_input("ecutwfc (Ry)", value=40.0, min_value=10.0, step=5.0)
with col3:
    ecutrho = st.number_input("ecutrho (Ry)", value=320.0, min_value=40.0, step=10.0)

kpoints_str = st.text_input("K-points grid (e.g., 4 4 4)", "4 4 4")
kpoints = kpoints_str.strip().split()

# Generate QE Input File
if st.button("🛠 Generate QE Input File"):
    if not atoms:
        st.error("Please upload a valid CIF file first.")
    elif len(kpoints) != 3 or not all(k.isdigit() for k in kpoints):
        st.error("Invalid k-point grid. Please provide 3 integers, e.g., 4 4 4")
    elif len(pseudo_map) < len(species):
        missing = [el for el in species if el not in pseudo_map]
        st.error(f"Missing pseudopotentials for: {', '.join(missing)}")
    else:
        qe_input, run_script = generate_qe_input(
            atoms, pseudo_map, calc_type, ecutwfc, ecutrho, kpoints
        )
        if qe_input:
            st.success("QE input file generated!")

            st.subheader("📄 QE Input Preview")
            st.code(qe_input, language="text")

            st.download_button("⬇ Download QE Input File", qe_input, file_name="qe_input.in")

            st.subheader("📜 QE Run Script (`pw.sh`)")
            st.code(run_script, language="bash")
            st.download_button("⬇ Download Run Script", run_script, file_name="pw.sh")
