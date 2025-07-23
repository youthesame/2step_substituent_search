from pathlib import Path

import pandas as pd
from ccdc import io
from ccdc.molecule import Atom, Bond
from tqdm import tqdm


def is_alkyl_chain_bond(bond, mol):
    """
    結合がアルキル鎖の一部かどうかを判定する

    Parameters
    ----------
    bond : Bond
        判定する結合
    mol : Molecule
        分子オブジェクト

    Returns
    -------
    bool
        アルキル鎖の一部の場合はTrue
    """
    a0, a1 = bond.atoms

    # C-C結合でない場合はFalse
    if a0.atomic_number != 6 or a1.atomic_number != 6:
        return False

    # 両方の炭素が4つの結合を持つ場合（sp3炭素）
    if len(a0.bonds) == 4 and len(a1.bonds) == 4:
        # 水素原子の数を数える
        h_count_a0 = sum(1 for neighbour in a0.neighbours if neighbour.atomic_number == 1)
        h_count_a1 = sum(1 for neighbour in a1.neighbours if neighbour.atomic_number == 1)

        # 両方の炭素が2つ以上の水素を持つ場合（-CH2-CH2-）
        if h_count_a0 >= 2 and h_count_a1 >= 2:
            return True

        # 一方の炭素が3つの水素を持つ場合（CH3-CH2-）
        if (h_count_a0 == 3 and h_count_a1 >= 2) or (h_count_a1 == 3 and h_count_a0 >= 2):
            return True

    return False


def fragment_molecule(molfile: str, max_ring: int = 8, exclude_alkyl_chains: bool = True) -> list:
    """
    指定された分子ファイルから分子を読み込み、特定の条件に基づいて分子をフラグメント化する。
    アルキル鎖除外ロジックを修正した版。

    Parameters
    ----------
    molfile : str
        読み込む分子のファイルパス。
    max_ring : int, optional
        保持する環の最大サイズ。デフォルトは8。
    exclude_alkyl_chains : bool, optional
        フラグメント化する際に、アルキル鎖を除外するかどうか。デフォルトはTrue。

    Returns
    -------
    list
        フラグメント化された分子のリスト。

    Raises
    ------
    RuntimeError
        分子の読み込みに失敗した場合に発生。
    """
    try:
        mol = io.MoleculeReader(molfile)
        mol = mol[0]
    except RuntimeError:
        print(f"RuntimeError: {molfile}")
        return None

    # 標準化（スクリーニング時に標準化しているなら、コメントアウトしてください）
    # mol.assign_bond_types()
    # mol.add_hydrogens()
    # mol.normalise_labels()

    # 切る予定の結合の集合
    weak_bonds = set()

    # 単結合を weak_bonds に追加
    for b in mol.bonds:
        if b.bond_type == 1:
            weak_bonds.add(b)
            # print(weak_bonds)

    # max_ring 以下の環に含まれる結合は weak_bonds から除外
    for r in mol.rings:
        if len(r.bonds) <= max_ring:
            for b in r.bonds:
                if b in weak_bonds:
                    weak_bonds.remove(b)

    removed_bonds = set()
    for b in weak_bonds:
        a0 = b.atoms[0]
        a1 = b.atoms[1]
        n0 = a0.atomic_number
        n1 = a1.atomic_number

        # 両方とも非炭素原子である場合は除外する。
        if n0 != 6 and n1 != 6:
            removed_bonds.add(b)

        # 末端水素は除外する。
        elif len(a0.bonds) == 1 and a0.atomic_symbol == "H":
            removed_bonds.add(b)
        elif len(a1.bonds) == 1 and a1.atomic_symbol == "H":
            removed_bonds.add(b)

        # 末端原子は除外する。
        elif len(a0.bonds) == 1 or len(a1.bonds) == 1:
            removed_bonds.add(b)

        # -CH2-CH2- は除外する。
        # elif n0 == 6 and n1 == 6:
        #     if len(a0.bonds) == 4 and len(a1.bonds) == 4:
        #         # 両方の炭素が2つの水素原子を持つ場合のみ除外
        #         h_count_a0 = sum(1 for neighbour in a0.neighbours if neighbour.atomic_number == 1)
        #         h_count_a1 = sum(1 for neighbour in a1.neighbours if neighbour.atomic_number == 1)
        #         if h_count_a0 == 2 and h_count_a1 == 2:
        #             removed_bonds.add(b)
        #             continue

        # アルキル鎖の除外
        elif exclude_alkyl_chains and is_alkyl_chain_bond(b, mol):
            removed_bonds.add(b)

    weak_bonds = weak_bonds - removed_bonds
    # print(weak_bonds)

    # 切った結合を保持しておくリスト
    atom_pairs = []

    # 結合を切る
    for b in weak_bonds:
        a0 = b.atoms[0]
        a1 = b.atoms[1]
        # 切る直前に atom_pairs に登録
        atom_pairs.append(b.atoms)

        # At(アスタチン、原子番号85)を付加する。
        try:
            at1 = Atom("At", coordinates=a0.coordinates)
            at2 = Atom("At", coordinates=a1.coordinates)
            at_0 = mol.add_atom(at1)
            at_1 = mol.add_atom(at2)
            mol.add_bond(Bond.BondType(1), b.atoms[0], at_1)
            mol.add_bond(Bond.BondType(1), b.atoms[1], at_0)
        except RuntimeError:
            print(mol.identifier)
        # 切る
        mol.remove_bond(b)

    return mol.components


def generate_smiles(mol: io.Molecule) -> tuple:
    """
    分子からSMILES表記を生成する。

    Parameters
    ----------
    mol : Molecule
        SMILES表記を生成する分子。

    Returns
    -------
    tuple
        AtのSMILES表記と、HのSMILES表記のタプル。
    """
    h_mol = mol.copy()
    h_mol = at_to_h(h_mol)
    return mol.smiles, h_mol.smiles


def at_to_h(mol: io.Molecule) -> io.Molecule:
    """
    分子内のアスタチン(At)原子を水素(H)原子に置換する。

    Parameters
    ----------
    mol : Molecule
        置換を行う分子。

    Returns
    -------
    Molecule
        アスタチンを水素に置換した後の分子。
    """
    # なぜか 1 回では At が消えないので、5 回繰り返す。
    for _ in range(5):
        for bond in mol.bonds:
            a0 = bond.atoms[0]
            a1 = bond.atoms[1]
            if a0.atomic_number == 85:
                h0 = Atom("H", coordinates=a0.coordinates)
                mol.add_bond(Bond.BondType(1), a1, (mol.add_atom(h0)))
                mol.remove_atom(a0)

            elif a1.atomic_number == 85:
                h1 = Atom("H", coordinates=a1.coordinates)
                mol.add_bond(Bond.BondType(1), a0, (mol.add_atom(h1)))
                mol.remove_atom(a1)

        for atom in mol.atoms:
            if atom.atomic_number == 85:
                continue
        break
    return mol


def include_metal(mol: io.Molecule) -> bool:
    """
    分子にAtを除いた金属原子が含まれているかどうかを判定する。

    Parameters
    ----------
    mol : Molecule
        判定する分子。

    Returns
    -------
    bool
        金属原子を含んでいる場合はTrue、含んでいない場合はFalse。
    """
    # kr(原子番号:36) 以降は DFT で計算できないため除外
    for atom in mol.atoms:
        if atom.atomic_number in (21, 22, 23, 24, 25, 26, 27, 28, 29, 30) or (
            atom.atomic_number >= 36 and atom.atomic_number != 85
        ):
            return True
    return False


def include_at(mol: io.Molecule) -> bool:
    """
    分子にAt原子が含まれているかどうかを判定する。

    Parameters
    ----------
    mol : Molecule
        判定する分子。

    Returns
    -------
    bool
        At原子を含んでいる場合はTrue、含んでいない場合はFalse。
    """
    for atom in mol.atoms:
        if atom.atomic_number == 85:
            return True
    return False


def include_nitoro(mol: io.Molecule) -> bool:
    """
    分子にニトロ基が含まれているかどうかを判定する。

    Parameters
    ----------
    mol : Molecule
        判定する分子。

    Returns
    -------
    bool
        ニトロ基を含んでいる場合はTrue、含んでいない場合はFalse。
    """
    atms = [a for a in mol.atoms if a.atomic_number == 7]
    # N がなければ、False を返す
    if len(atms) == 0:
        return False
    else:
        for a in atms:
            bnds = a.bonds
            n2 = []
            counts = 0
            for b in bnds:
                if b.bond_type == 2:
                    n2.append(b)
            # N が二重結合を2つ持ち、両方とも結合している原子が O であれば True を返す
            if len(n2) == 2:
                for n in n2:
                    a0 = n.atoms[0]
                    a1 = n.atoms[1]
                    if a0.atomic_number == 8 or a1.atomic_number == 8:
                        counts += 1
                if counts == 2:
                    return True
        return False


def test_fragmentation():
    molfile = "../../data/test/fragmentation/POLFAB.mol"
    fragments = fragment_molecule(molfile, exclude_alkyl_chains=True)
    # 出力ディレクトリの作成
    output_at_dir = Path("../../data/test/fragment_mol_at")
    output_h_dir = Path("../../data/test/fragment_mol_h")
    output_at_dir.mkdir(parents=True, exist_ok=True)
    output_h_dir.mkdir(parents=True, exist_ok=True)

    existing_at_smiles = {}
    existing_h_smiles = {}
    f_count = 0
    h_count = 0

    for fragment in fragments:
        at_smiles, h_smiles = generate_smiles(fragment)
        if at_smiles not in existing_at_smiles:
            existing_at_smiles[at_smiles] = f_count
            f_count += 1
            at_output_path = output_at_dir / f"{f_count}_fragment.mol"
            with io.MoleculeWriter(str(at_output_path)) as f:
                f.write(fragment)
        if h_smiles not in existing_h_smiles:
            existing_h_smiles[h_smiles] = h_count
            h_count += 1
            h_output_path = output_h_dir / f"{h_count}_h_fragment.mol"
            with io.MoleculeWriter(str(h_output_path)) as f:
                f.write(fragment)


def main():
    # 出力ディレクトリの作成
    output_at_dir = Path("../../data/substituent/fragment_mol_at")
    output_h_dir = Path("../../data/substituent/fragment_mol_h")
    output_at_dir.mkdir(parents=True, exist_ok=True)
    output_h_dir.mkdir(parents=True, exist_ok=True)

    molfiles = sorted(Path("../../data/substituent/mol_original").glob("*_original.mol"))
    fragments_data = []
    existing_at_smiles = {}
    existing_h_smiles = {}
    f_count = 0
    h_count = 0
    # glob した分子ファイルを順番に読み込む
    for molfile in tqdm(molfiles, total=len(molfiles)):
        try:
            fragments = fragment_molecule(str(molfile))
            if fragments is None:
                continue

            for fragment in fragments:
                # フラグメント化できなかった分子、つまり At を持たない分子は除外
                if include_at(fragment):
                    # At を除く金属原子を含む分子は除外
                    if not include_metal(fragment):
                        at_smiles, h_smiles = generate_smiles(fragment)

                        # At の SMILES が重複している分子は除外
                        if at_smiles not in existing_at_smiles:
                            f_code = f"{f_count:06d}"
                            existing_at_smiles[at_smiles] = f_code
                            f_count += 1

                            # At 付き分子を保存
                            at_output_path = output_at_dir / f"{f_code}_fragment.mol"
                            with io.MoleculeWriter(str(at_output_path)) as f:
                                f.write(fragment)

                            # At を H に置き換えた分子を保存、重複の場合は保存しない
                            if h_smiles not in existing_h_smiles:
                                h_fragment = at_to_h(fragment)
                                h_fragment.normalise_hydrogens()
                                h_code = f"{h_count:06d}"
                                existing_h_smiles[h_smiles] = h_code
                                h_count += 1
                                h_output_path = output_h_dir / f"{h_code}_h_fragment.mol"
                                with io.MoleculeWriter(str(h_output_path)) as f:
                                    f.write(h_fragment)

                            else:
                                h_code = existing_h_smiles[h_smiles]

                            # 分子情報の保存
                            fragments_data.append([
                                existing_at_smiles[at_smiles],
                                h_code,
                                at_smiles,
                                h_smiles,
                                len(h_fragment.atoms),
                                h_fragment.molecular_weight,
                            ])

        except Exception as e:
            print(f"Error: {molfile}")
            print(e)
            continue

    df = pd.DataFrame(
        fragments_data,
        columns=[
            "fragment_code",
            "h_code",
            "At_SMILES",
            "H_SMILES",
            "n_atoms",
            "weight",
        ],
    )
    df["n_At"] = df["At_SMILES"].str.count("At")
    df.to_csv("../../data/substituent/fragments_data.csv", index=False)


# TODO: H_SMILESが重複している分子を取り出して、１つの分子に At を集約する？
# At の個数などを気にする必要あり。

if __name__ == "__main__":
    main()
    # test_fragmentation()
    # fragment_codeとh_codeはstr型に変換して"fragments_data.csv"読み込む
    # df = pd.read_csv("fragments_data.csv", dtype={"fragment_code": str, "h_code": str})
    # print(df)
