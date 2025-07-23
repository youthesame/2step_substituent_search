import time
from datetime import datetime
from pathlib import Path

import pandas as pd
from ccdc import CSD_VERSION_LATEST, io
from ccdc.search import Search
from tqdm import tqdm

"""
screen_single_molecule_crystals()
Screening for single-molecule crystals:
- Only crystals containing a single unique molecule (len(set(SMILES)) == 1)
- Has 3D coordinates
- No errors
- Not from powder X-ray diffraction
- Not polymeric
- Records REFCODE and unit cell parameters (a, b, c, α, β, γ)
"""

# Define the base directory for data/substituent at the project root
BASE_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "substituent"


def record_csd_version(md_file_path: str = "screening_info.md", elapsed_time: float = None):
    """
    Records the current date, CSD version, and execution time in a Markdown file.

    Parameters
    ----------
    md_file_path : str
        Path to the Markdown file to record information.
    elapsed_time : float, optional
        Execution time in seconds.
    """
    current_date = datetime.now().strftime("%Y-%m-%d")
    csd_version = CSD_VERSION_LATEST

    content = f"# Screening Information\n\n- Execution Date: {current_date}\n- CSD version used: {csd_version}\n"
    if elapsed_time is not None:
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = int(elapsed_time % 60)
        content += f"- Execution Time: {hours}h {minutes}m {seconds}s\n"

    with Path(md_file_path).open("w", encoding="utf-8", newline="\n") as md_file:
        md_file.write(content)


def screen_single_molecule_crystals():
    """
    Screens for single-molecule crystals from the CSD and records unit cell parameters.

    Screening criteria:
    - Only crystals with a single unique molecule (len(set(SMILES)) == 1)
    - Has 3D coordinates
    - No errors
    - Not from powder X-ray diffraction
    - Not polymeric

    Output: parquet file with SMILES, InChI, unit cell parameters (a,b,c,α,β,γ), and REFCODE
    """
    record_csd_version()

    # Set up screening conditions
    setting = Search.Settings()
    setting.has_3d_coordinates = True
    setting.no_errors = True
    setting.no_powder = True
    setting.not_polymeric = True

    crystal_data = []

    with io.EntryReader("CSD") as reader:
        for entry in tqdm(reader, desc="Screening CSD entries"):
            try:
                # Skip entries that do not meet the basic settings criteria
                if not setting.test(entry):
                    continue

                # Get all molecule components and their SMILES
                smiles_list = []
                molecules = []

                for mol in entry.molecule.components:
                    # Complete missing bonds and hydrogens
                    mol.assign_bond_types(which="unknown")
                    mol.add_hydrogens(mode="missing")

                    # Get SMILES; skip if error
                    smiles = mol.smiles
                    if smiles is None:
                        continue

                    smiles_list.append(smiles)
                    molecules.append(mol)

                # Check if this is a single-molecule crystal
                unique_smiles = set(smiles_list)
                if len(unique_smiles) != 1:
                    continue  # Skip multi-molecule crystals

                # Get the single unique molecule
                representative_mol = molecules[0]
                unique_smiles_str = list(unique_smiles)[0]

                # Generate InChI
                try:
                    inchi = representative_mol.generate_inchi().inchi
                except Exception:
                    inchi = None

                # Get unit cell parameters
                crystal = entry.crystal
                if crystal is None:
                    continue

                cell_lengths = crystal.cell_lengths  # (a, b, c)
                cell_angles = crystal.cell_angles  # (α, β, γ)

                # Record the data
                crystal_data.append({
                    "SMILES": unique_smiles_str,
                    "InChI": inchi,
                    "a": cell_lengths.a,
                    "b": cell_lengths.b,
                    "c": cell_lengths.c,
                    "alpha": cell_angles.alpha,
                    "beta": cell_angles.beta,
                    "gamma": cell_angles.gamma,
                    "refcode": entry.identifier,
                })

            except Exception as e:
                print(f"Error processing {entry.identifier}: {e}")
                continue

    # Create DataFrame
    df = pd.DataFrame(crystal_data)
    BASE_DIR.mkdir(parents=True, exist_ok=True)

    # Save as Parquet
    parquet_file = BASE_DIR / "screening_refcode.parquet"
    df.to_parquet(parquet_file, index=False)

    print(f"Screening complete. Found {len(df)} single-molecule crystals.")

    return df


if __name__ == "__main__":
    start_time = time.time()
    screen_single_molecule_crystals()
    elapsed_time = time.time() - start_time
    record_csd_version(elapsed_time=elapsed_time)
