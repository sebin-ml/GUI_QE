import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from ase.io import read
from ase.data import atomic_masses, chemical_symbols
import os
import shutil

# Helper: get atomic number from symbol
def ase_symbol_to_z(symbol):
    return chemical_symbols.index(symbol)

class CIF2QEConverter:
    def __init__(self, root):
        self.root = root
        self.root.title("CIF to QE Input Converter")
        self.root.geometry("600x400")
        self.pseudopotentials = {}
        self.atoms = None
        self.species = []

        # CIF File Selection
        tk.Label(root, text="Step 1: Select CIF file", font=("Arial", 12)).pack(pady=10)
        tk.Button(root, text="Browse CIF", command=self.select_cif, font=("Arial", 11), bg="lightblue").pack()
        self.status_label = tk.Label(root, text="", font=("Arial", 10))
        self.status_label.pack(pady=5)

        # QE Parameters
        frame = tk.Frame(root)
        frame.pack(pady=10)

        tk.Label(frame, text="Calculation:", font=("Arial", 10)).grid(row=0, column=0, padx=5, sticky="e")
        self.calc_type = tk.Entry(frame, width=10)
        self.calc_type.insert(0, "scf")
        self.calc_type.grid(row=0, column=1)

        tk.Label(frame, text="ecutwfc:", font=("Arial", 10)).grid(row=0, column=2, padx=5, sticky="e")
        self.ecutwfc = tk.Entry(frame, width=10)
        self.ecutwfc.insert(0, "40")
        self.ecutwfc.grid(row=0, column=3)

        tk.Label(frame, text="ecutrho:", font=("Arial", 10)).grid(row=0, column=4, padx=5, sticky="e")
        self.ecutrho = tk.Entry(frame, width=10)
        self.ecutrho.insert(0, "320")
        self.ecutrho.grid(row=0, column=5)

        tk.Label(frame, text="K-Points (e.g., 4 4 4):", font=("Arial", 10)).grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="e")
        self.kpoints = tk.Entry(frame, width=20)
        self.kpoints.insert(0, "4 4 4")
        self.kpoints.grid(row=1, column=2, columnspan=3)

        # Proceed button
        tk.Button(root, text="Continue to Pseudopotential Selection", command=self.ask_pseudopotentials,
                  bg="lightgreen", font=("Arial", 12)).pack(pady=20)

    def select_cif(self):
        cif_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif")])
        if not cif_path:
            return

        try:
            self.atoms = read(cif_path)
            self.species = sorted(set(self.atoms.get_chemical_symbols()))
            self.cif_path = cif_path
            self.status_label.config(text=f"Detected elements: {', '.join(self.species)}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read CIF file:\n{str(e)}")

    def ask_pseudopotentials(self):
        if not self.atoms:
            messagebox.showerror("Error", "Please select a CIF file first.")
            return

        self.popup = tk.Toplevel(self.root)
        self.popup.title("Select Pseudopotentials")
        self.popup.geometry("600x100")

        self.labels = {}

        for i, element in enumerate(self.species):
            tk.Label(self.popup, text=f"{element}:", font=("Arial", 10)).grid(row=i, column=0, padx=10, pady=5, sticky="e")
            label = tk.Label(self.popup, text="No file selected", width=40, anchor="w")
            label.grid(row=i, column=1, padx=10)
            self.labels[element] = label

            tk.Button(self.popup, text="Select .UPF", command=lambda e=element: self.select_upf(e)).grid(row=i, column=2)

        tk.Button(self.popup, text="Generate QE Input File", command=self.generate_input, bg="lightgreen", font=("Arial", 11)).grid(
            row=len(self.species), column=1, pady=20)

    def select_upf(self, element):
        upf_path = filedialog.askopenfilename(filetypes=[("UPF files", "*.UPF")])
        if upf_path:
            filename = os.path.basename(upf_path)
            self.pseudopotentials[element] = filename
            self.labels[element].config(text=filename)

            if not os.path.exists(filename):
                try:
                    shutil.copy(upf_path, filename)
                except Exception as e:
                    messagebox.showwarning("Copy Failed", f"Could not copy {filename}:\n{str(e)}")

    def generate_input(self):
        missing = [el for el in self.species if el not in self.pseudopotentials]
        if missing:
            messagebox.showerror("Missing Pseudopotentials", f"Please select UPF files for: {', '.join(missing)}")
            return

        out_file = os.path.splitext(self.cif_path)[0] + ".in"
        run_file = os.path.splitext(self.cif_path)[0] + "_pw.sh"

        try:
            lines = []
            lines.append("&control")
            lines.append(f"  calculation = '{self.calc_type.get()}',")
            lines.append("  prefix = 'calc',")
            lines.append("  outdir = './tmp',")
            lines.append("/\n")

            lines.append("&system")
            lines.append(f"  ibrav = 0,")
            lines.append(f"  nat = {len(self.atoms)},")
            lines.append(f"  ntyp = {len(set(self.atoms.get_chemical_symbols()))},")
            lines.append(f"  ecutwfc = {self.ecutwfc.get()},")
            lines.append(f"  ecutrho = {self.ecutrho.get()},")
            lines.append("/\n")

            lines.append("&electrons")
            lines.append("  conv_thr = 1.0d-6")
            lines.append("/\n")

            lines.append("CELL_PARAMETERS angstrom")
            for vec in self.atoms.get_cell():
                lines.append("  {:.10f} {:.10f} {:.10f}".format(*vec))

            lines.append("\nATOMIC_SPECIES")
            for el in self.species:
                z = ase_symbol_to_z(el)
                mass = atomic_masses[z]
                lines.append(f"{el}  {mass:.4f}  {self.pseudopotentials[el]}")

            lines.append("\nATOMIC_POSITIONS angstrom")
            for s in self.atoms:
                lines.append(f"{s.symbol}  {s.position[0]:.6f}  {s.position[1]:.6f}  {s.position[2]:.6f}")

            lines.append("\nK_POINTS automatic")
            kpts = self.kpoints.get().strip().split()
            if len(kpts) != 3:
                messagebox.showerror("Invalid K-Points", "Please enter 3 integers for k-points (e.g., 4 4 4)")
                return
            lines.append("  " + " ".join(kpts) + "  0 0 0\n")

            with open(out_file, 'w') as f:
                f.write("\n".join(lines))

            # Write pw.sh script
            with open(run_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("# Run Quantum ESPRESSO pw.x calculation\n")
                f.write(f"pw.x < {os.path.basename(out_file)} > output.log\n")

            os.chmod(run_file, 0o755)  # make executable

            self.preview_file(out_file)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to write files:\n{str(e)}")

    def preview_file(self, filepath):
        try:
            with open(filepath, 'r') as f:
                content = f.read()

            preview = tk.Toplevel(self.root)
            preview.title("Preview of QE Input File")
            preview.geometry("700x600")

            text = scrolledtext.ScrolledText(preview, wrap=tk.WORD, font=("Courier", 10))
            text.insert(tk.END, content)
            text.pack(expand=True, fill='both')

        except Exception as e:
            messagebox.showerror("Preview Error", str(e))


# Run GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = CIF2QEConverter(root)
    root.mainloop()
