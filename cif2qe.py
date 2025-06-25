import streamlit as st
import os
import tempfile
from ase.io import read
from ase.data import atomic_masses, chemical_symbols
import shutil
import io

def ase_symbol_to_z(symbol):
    return chemical_symbols.index(symbol)

def generate_qe_input(cif_file, pseudo_files, calc_type, ecutwfc, ecutrho, kpoints):
    atoms = read(cif_file)
    species = sorted(set(atoms.get_chemical_symbols()))

    # Save pseudopotentials temporarily
    with tempfile.TemporaryDirectory() as tmpdir:
        pseudo_dir = os.path.join(tmpdir, "pseudos")
        os.makedirs(pseudo_dir, exist_ok=True)

        pseudo_map = {}
        for file in pseudo_files:
            for el in species:
                if file.name.lower().startswith(el.lower()) and file.name.lower().endswith('.upf'):
                    filepath = os.path.join(pseudo_dir, file.name)
                    with open(filepath, "wb") as f:
                        f.write(file.read())
                    pseudo_map[el] = file.name

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
            lines.append(f"{el}  {mass:.4f}  {pseudo_map[el]}")

        lines.append("\nATOMIC_POSITIONS angstrom")
        for atom in atoms:
            pos = atom.position
            lines.append(f"{atom.symbol}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")

        lines.append("\nK_POINTS automatic")
        lines.append("  " + " ".join(kpoints) + " 0 0 0")

        # Create input file and run script
        input_str = "\n".join(lines)
        run_script = f"""#!/bin/bash
# Run Quantum ESPRESSO
pw.x < qe_input.in > output.log
"""

        return input_str, run_script


# ---------- STREAMLIT UI ----------
st.title("🧪 CIF to QE Input File Converter")
st.markdown("Upload your CIF file, select pseudopotentials, and configure Quantum ESPRESSO parameters.")

# Step 1: Upload CIF
cif_file = st.file_uploader("📁 Upload CIF file", type=["cif"])

# Step 2: Upload pseudopotentials
pseudo_files = st.file_uploader("📤 Upload .UPF pseudopotentials", type=["upf", "UPF"], accept_multiple_files=True)

# Step 3: QE parameters
col1, col2, col3 = st.columns(3)
with col1:
    calc_type = st.selectbox("Calculation Type", ["scf", "relax", "vc-relax", "nscf"])
with col2:
    ecutwfc = st.number_input("ecutwfc (Ry)", min_value=10.0, value=40.0)
with col3:
    ecutrho = st.number_input("ecutrho (Ry)", min_value=40.0, value=320.0)

kpoints = st.text_input("K-Points Grid (e.g. 4 4 4)", "4 4 4")

# Step 4: Generate input file
if st.button("⚙ Generate QE Input"):
    if not cif_file or not pseudo_files:
        st.warning("Please upload both CIF and pseudopotential files.")
    else:
        kpt_list = kpoints.strip().split()
        if len(kpt_list) != 3 or not all(k.isdigit() for k in kpt_list):
            st.error("❌ Invalid k-points. Enter 3 integers like: 4 4 4")
        else:
            qe_input, run_script = generate_qe_input(
                cif_file,
                pseudo_files,
                calc_type,
                ecutwfc,
                ecutrho,
                kpt_list
            )

            if qe_input:
                st.success("✅ QE input file generated!")

                st.subheader("📄 Preview: QE Input File")
                st.code(qe_input, language="bash")

                st.download_button("⬇ Download QE Input File", qe_input, file_name="qe_input.in")

                st.subheader("🛠 Run Script (pw.sh)")
                st.code(run_script, language="bash")
                st.download_button("⬇ Download Run Script", run_script, file_name="pw.sh")

