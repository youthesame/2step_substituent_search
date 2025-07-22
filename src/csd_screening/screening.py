import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
from ccdc import CSD_VERSION_LATEST, io
from ccdc.search import Search
from tqdm import tqdm

"""
prescreening_molecules()
スクリーニング step1
- 元素は6-31Gを適用可能なH-Krまで
- 3D座標を持つ
- errorのない
- イオン性のものは除く
- 粉末X線回折で構造解析されているものは除く
- 高分子は除く
- ラジカルを含まない

filter_and_reassign_molcodes()
スクリーニング step2
- SMILES 重複したら、R値が小さいものを優先する。
"""


def record_csd_version(md_file_path: str = "screening_info.md"):
    """
    現在の日付とCSDのバージョンをMarkdownファイルに記録します。

    Parameters
    ----------
    md_file_path : str
        情報を記録するMarkdownファイルのパス。
    """
    current_date = datetime.now().strftime("%Y-%m-%d")
    csd_version = CSD_VERSION_LATEST

    content = f"# Screening Information\n\n- Execution Date: {current_date}\n- CSD version used: {csd_version}\n"

    with Path(md_file_path).open("w", encoding="utf-8") as md_file:
        md_file.write(content)


def prescreening_molecules():
    """
    化学データベースから分子をスクリーニングし、MOLファイルと基本情報を保存します。

    以下の条件を満たす分子をスクリーニングします：
    - 元素は6-31Gを適用可能なH-Krまで
    - 3D座標を持つ
    - エラーがない
    - イオン性のものは除く
    - 粉末X線回折で構造解析されているものは除く
    - 高分子は除く
    - ラジカルを含まない
    """
    record_csd_version()
    setting = Search.Settings()
    setting.has_3d_coordinates = True
    setting.no_errors = True
    setting.no_ions = True
    setting.no_powder = True
    setting.not_polymeric = True
    # setting.no_disorder = True
    # setting.no_metals = True
    # setting.only_organic = True
    # setting.only_organometallic = False

    mol_data = []
    all_smiles = []
    r_factors = []
    with io.EntryReader("CSD") as reader:
        for entry in tqdm(reader):
            try:
                # settingの条件を満たさないものはスキップ
                if not setting.test(entry):
                    continue

                r_factor = entry.r_factor
                disorder = entry.has_disorder
                for mol in entry.molecule.components:
                    # 欠損の補完
                    mol.assign_bond_types(which="unknown")
                    mol.add_hydrogens(mode="missing")

                    # SMILESの取得、エラーならスキップ
                    smiles = mol.smiles
                    if smiles is None:
                        continue

                    # kr(原子番号:36) 以降はスキップ
                    if any(atom.atomic_number >= 36 for atom in mol.atoms):
                        continue

                    # ラジカルを含むものはスキップ
                    electrons = sum([a.atomic_number for a in mol.atoms])
                    if electrons % 2 == 1:
                        continue

                    try:
                        molcode = all_smiles.index(smiles)
                        new = False
                    except ValueError:
                        molcode = len(all_smiles)
                        new = True

                    # MOLファイルの保存
                    if (
                        new
                        or r_factors[molcode] is None
                        or (r_factor is not None and r_factor < r_factors[molcode])
                    ):
                        mol_path = Path(f"prescreening/{molcode:06d}_prescreening.mol")
                        mol_path.parent.mkdir(parents=True, exist_ok=True)
                        with mol_path.open("w") as mol_file:
                            mol_file.write(mol.to_string(format="mol"))

                        if new:
                            all_smiles.append(smiles)
                            r_factors.append(r_factor)
                        else:
                            r_factors[molcode] = r_factor

                        # データの保存
                        molcode = f"{molcode:06d}"
                        mol_data.append(
                            (
                                molcode,
                                smiles,
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

    # parquetファイルに保存する
    df = pd.DataFrame(
        mol_data,
        columns=[
            "molcode",
            "SMILES",
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
    df.to_parquet("screening_step1.parquet")


def filter_and_update_molcodes(parquet_file: str, output_csv: str):
    """
    スクリーニングされた分子のデータから、SMILES表記が重複している分子の中でR値が最小のものを選択し、molcodeを更新して保存する。

    Parameters
    ----------
    parquet_file : str
        入力parquetファイルのパス
    output_csv : str
        出力CSVファイルのパス
    """
    df = pd.read_parquet(parquet_file)

    # r_factor列でNaN値を大きな数値に置き換える
    df["r_factor"].fillna(1e10, inplace=True)

    # SMILESが重複している分子の中でR値が最小のものを選択
    filtered_df = df.loc[df.groupby("SMILES")["r_factor"].idxmin()]

    Path("mol_original").mkdir(parents=True, exist_ok=True)
    # molcodeを更新
    for idx, row in tqdm(enumerate(filtered_df.itertuples()), total=len(filtered_df)):
        new_molcode = f"{idx:06d}"
        old_mol_path = Path(f"prescreening/{row.molcode}_prescreening.mol")
        new_mol_path = Path(f"mol_original/{new_molcode}_original.mol")
        shutil.move(old_mol_path, new_mol_path)

        filtered_df.at[row.Index, "molcode"] = new_molcode

    # csvファイルに保存する
    filtered_df.to_csv(output_csv, index=False)


def archive_del_dir(directory_path: str):
    """
    指定されたディレクトリをZIPファイルにアーカイブし、その後ディレクトリを削除します。

    Parameters
    ----------
    directory_path : str
        アーカイブするディレクトリのパス。
    """
    directory = Path(directory_path)

    # ディレクトリが存在しない場合は何もしない
    if not directory.is_dir():
        print(f"ディレクトリが存在しません: {directory_path}")
        return

    # ディレクトリをZIPファイルにアーカイブ
    shutil.make_archive(directory_path, "zip", directory.parent, directory.name)

    # 元のディレクトリを削除
    shutil.rmtree(directory)


if __name__ == "__main__":
    prescreening_molecules()
    filter_and_update_molcodes("screening_step1.parquet", "screening_step2.csv")
    archive_del_dir("prescreening")
