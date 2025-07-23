import shutil
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
from ccdc import CSD_VERSION_LATEST, io
from ccdc.search import Search
from tqdm import tqdm

"""
prescreening_molecules()
Screening step 1
- Elements up to H-Kr (suitable for 6-31G)
- Has 3D coordinates
- No errors
- Excludes ionic compounds
- Excludes structures determined by powder X-ray diffraction
- Excludes polymers
- Excludes radicals

filter_and_reassign_molcodes()
Screening step 2
- If SMILES are duplicated, prioritize the one with the smallest R value.
"""

# Define the base directory for data/substituent at the project root
BASE_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "substituent"


def record_csd_version(
    md_file_path: str = "screening_info.md", elapsed_time: float = None
):
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


def prescreening_molecules():
    """
    Screens molecules from a chemical database and saves MOL files and basic information.

    The following criteria are used for screening:
    - Elements up to H-Kr (suitable for 6-31G)
    - Has 3D coordinates
    - No errors
    - Excludes ionic compounds
    - Excludes structures determined by powder X-ray diffraction
    - Excludes polymers
    - Excludes radicals
    """
    record_csd_version()
    setting = Search.Settings()
    setting.has_3d_coordinates = True
    setting.no_errors = True
    setting.no_ions = True
    setting.no_powder = True
    setting.not_polymeric = True
    setting.no_metals = True
    setting.only_organic = True
    # setting.only_organometallic = False

    mol_data = []
    all_smiles = []
    r_factors = []
    with io.EntryReader("CSD") as reader:
        for entry in tqdm(reader):
            try:
                # Skip entries that do not meet the settings criteria
                if not setting.test(entry):
                    continue

                r_factor = entry.r_factor
                disorder = entry.has_disorder
                for mol in entry.molecule.components:
                    # Complete missing bonds and hydrogens
                    mol.assign_bond_types(which="unknown")
                    mol.add_hydrogens(mode="missing")

                    # Get SMILES; skip if error
                    smiles = mol.smiles
                    if smiles is None:
                        continue

                    # Skip if any atom is Kr (atomic number >= 36) or heavier
                    if any(atom.atomic_number >= 36 for atom in mol.atoms):
                        continue

                    # Skip if contains radicals (odd number of electrons)
                    electrons = sum([a.atomic_number for a in mol.atoms])
                    if electrons % 2 == 1:
                        continue

                    # Generate InChI
                    try:
                        inchi = mol.generate_inchi().inchi
                    except Exception:
                        inchi = None

                    try:
                        molcode = all_smiles.index(smiles)
                        new = False
                    except ValueError:
                        molcode = len(all_smiles)
                        new = True

                    # Save MOL file
                    if (
                        new
                        or r_factors[molcode] is None
                        or (r_factor is not None and r_factor < r_factors[molcode])
                    ):
                        mol_path = (
                            BASE_DIR
                            / "prescreening"
                            / f"{molcode:06d}_prescreening.mol"
                        )
                        mol_path.parent.mkdir(parents=True, exist_ok=True)
                        with mol_path.open("w") as mol_file:
                            mol_file.write(mol.to_string(format="mol"))

                        if new:
                            all_smiles.append(smiles)
                            r_factors.append(r_factor)
                        else:
                            r_factors[molcode] = r_factor

                        # Save data
                        molcode = f"{molcode:06d}"
                        mol_data.append(
                            (
                                molcode,
                                smiles,
                                inchi,
                                entry.identifier,
                                entry.chemical_name,
                                mol.formula,
                                mol.molecular_weight,
                                mol.molecular_volume,
                                electrons,
                                r_factor,
                                disorder,
                            )
                        )

            except Exception as e:
                print(f"Error: {entry.identifier}")
                print(e)
                continue

    df = pd.DataFrame(
        mol_data,
        columns=[
            "molcode",
            "SMILES",
            "InChI",
            "refcode",
            "entry_name",
            "formula",
            "weight",
            "volume",
            "electrons",
            "r_factor",
            "disorder",
        ],
    )
    (BASE_DIR).mkdir(parents=True, exist_ok=True)
    df.to_csv(BASE_DIR / "screening_step1.csv")


def filter_and_update_molcodes(csv_file: str, output_csv: str):
    """
    From the screened molecule data, select the molecule with the smallest R value among those with duplicate SMILES,
    update the molcode, and save.

    Parameters
    ----------
    csv_file : str
        Path to the input CSV file.
    output_csv : str
        Path to the output CSV file.
    """
    df = pd.read_csv(csv_file)
    df["molcode"] = df["molcode"].astype(str)
    df["r_factor"] = df["r_factor"].fillna(1e10)
    filtered_df = df.loc[df.groupby("SMILES")["r_factor"].idxmin()]
    filtered_df = filtered_df.reset_index(drop=True)

    (BASE_DIR / "mol_original").mkdir(parents=True, exist_ok=True)

    # Extract only rows for which the corresponding MOL file exists
    exists_mask = filtered_df["molcode"].apply(
        lambda mc: (BASE_DIR / "prescreening" / f"{mc}_prescreening.mol").exists()
    )
    missing_count = (~exists_mask).sum()
    filtered_df = filtered_df[exists_mask].reset_index(drop=True)

    # Reassign new molcodes and copy MOL files
    for idx, row in filtered_df.iterrows():
        new_molcode = f"{idx:06d}"
        old_mol_path = BASE_DIR / "prescreening" / f"{row['molcode']}_prescreening.mol"
        new_mol_path = BASE_DIR / "mol_original" / f"{new_molcode}_original.mol"
        shutil.copy(old_mol_path, new_mol_path)
        filtered_df.at[idx, "molcode"] = new_molcode

    if missing_count > 0:
        print(f"Total skipped rows due to missing MOL files: {missing_count}")

    filtered_df.to_csv(output_csv, index=False)


def archive_del_dir(directory_path: str):
    """
    Archives the specified directory as a ZIP file and then deletes the directory.

    Parameters
    ----------
    directory_path : str
        Path to the directory to archive.
    """
    directory = Path(directory_path)

    # Do nothing if the directory does not exist
    if not directory.is_dir():
        print(f"Directory does not exist: {directory_path}")
        return

    # Archive the directory as a ZIP file
    shutil.make_archive(directory_path, "zip", directory.parent, directory.name)

    # Delete the original directory
    shutil.rmtree(directory)


if __name__ == "__main__":
    start_time = time.time()
    prescreening_molecules()
    filter_and_update_molcodes(
        str(BASE_DIR / "screening_step1.csv"), str(BASE_DIR / "screening_step2.csv")
    )
    archive_del_dir(str(BASE_DIR / "prescreening"))
    elapsed_time = time.time() - start_time
    record_csd_version(elapsed_time=elapsed_time)
