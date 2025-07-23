from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ccdc import descriptors, io
from matplotlib.backends.backend_pdf import PdfPages


def analyze_fragments():
    fragment_dir = Path("../../data/substituent/fragment_mol_h")
    output_dir = Path("../../result/fragment")
    output_dir.mkdir(parents=True, exist_ok=True)

    mol_files = list(fragment_dir.glob("*.mol"))[:1000]

    atom_counts = []
    non_h_atom_counts = []
    atom_types = Counter()
    volumes = []

    for mol_file in mol_files:
        try:
            mol = io.MoleculeReader(str(mol_file))[0]

            atom_count = len(mol.atoms)
            atom_counts.append(atom_count)

            non_h_count = sum(1 for atom in mol.atoms if atom.atomic_symbol != "H")
            non_h_atom_counts.append(non_h_count)

            for atom in mol.atoms:
                atom_types[atom.atomic_symbol] += 1

            try:
                volume = descriptors.MolecularDescriptors.van_der_waals_volume(mol)
                if volume is not None:
                    volumes.append(volume)
            except:
                pass

        except Exception as e:
            print(f"Error processing {mol_file}: {e}")
            continue

    with PdfPages(output_dir / "fragment_analysis.pdf") as pdf:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

        ax1.hist(atom_counts, bins=30, alpha=0.7, edgecolor="black")
        ax1.set_xlabel("Total Atom Count")
        ax1.set_ylabel("Frequency")

        ax2.hist(non_h_atom_counts, bins=20, alpha=0.7, edgecolor="black")
        ax2.set_xlabel("Non-H Atom Count")
        ax2.set_ylabel("Frequency")

        common_atoms = atom_types.most_common(10)
        elements, counts = zip(*common_atoms) if common_atoms else ([], [])
        ax3.bar(elements, counts)
        ax3.set_xlabel("Atom Type")
        ax3.set_ylabel("Total Count")
        ax3.tick_params(axis="x", rotation=45)

        if volumes:
            ax4.hist(volumes, bins=25, alpha=0.7, edgecolor="black")
            ax4.set_xlabel("Van der Waals Volume (r)")
            ax4.set_ylabel("Frequency")
        else:
            ax4.text(0.5, 0.5, "No volume data available", transform=ax4.transAxes, ha="center", va="center")

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

    print(f"Analysis complete. Results saved to {output_dir / 'fragment_analysis.pdf'}")
    print(f"Processed {len(mol_files)} fragment files")
    print(f"Average atom count: {np.mean(atom_counts):.1f}")
    print(f"Average non-H atom count: {np.mean(non_h_atom_counts):.1f}")
    print(f"Most common atom types: {dict(atom_types.most_common(5))}")
    if volumes:
        print(f"Average volume: {np.mean(volumes):.1f} r")


if __name__ == "__main__":
    analyze_fragments()
