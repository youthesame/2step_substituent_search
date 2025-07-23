from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def find_matching_pairs(df: pd.DataFrame, axis1: str, axis2: str, tolerance: float = 0.25) -> pd.DataFrame:
    """
    Find pairs of crystals where two specified axes match within tolerance.

    Args:
        df: DataFrame containing crystal data
        axis1: First axis to compare ('a', 'b', or 'c')
        axis2: Second axis to compare ('a', 'b', or 'c')
        tolerance: Tolerance range for matching (±tolerance)

    Returns:
        DataFrame containing matching pairs with their differences
    """
    matches = []

    # Get unique combinations of crystal indices
    for i, j in combinations(range(len(df)), 2):
        row1 = df.iloc[i]
        row2 = df.iloc[j]

        # Check if both axes are within tolerance
        diff1 = abs(row1[axis1] - row2[axis1])
        diff2 = abs(row1[axis2] - row2[axis2])

        if diff1 <= tolerance and diff2 <= tolerance:
            matches.append({
                "refcode1": row1["refcode"],
                "refcode2": row2["refcode"],
                f"{axis1}_1": row1[axis1],
                f"{axis1}_2": row2[axis1],
                f"{axis2}_1": row1[axis2],
                f"{axis2}_2": row2[axis2],
                f"diff_{axis1}": diff1,
                f"diff_{axis2}": diff2,
                "a_1": row1["a"],
                "b_1": row1["b"],
                "c_1": row1["c"],
                "a_2": row2["a"],
                "b_2": row2["b"],
                "c_2": row2["c"],
            })

    return pd.DataFrame(matches)


def find_specific_value_matches(
    df: pd.DataFrame, target_values: List[Tuple[float, float]], tolerance: float = 0.25
) -> Dict[str, pd.DataFrame]:
    """
    Find crystals that match specific target value pairs within tolerance.

    Args:
        df: DataFrame containing crystal data
        target_values: List of (value1, value2) tuples to search for
        tolerance: Tolerance range for matching (±tolerance)

    Returns:
        Dictionary mapping target value strings to DataFrames of matches
    """
    results = {}

    for val1, val2 in target_values:
        print(f"Searching for matches with values {val1}, {val2}...")
        matches = []

        # Create boolean masks for efficient filtering
        axes = ["a", "b", "c"]

        # Check all combinations of axes
        for i in range(3):
            for j in range(i + 1, 3):
                axis1, axis2 = axes[i], axes[j]

                # Check both orientations: (val1, val2) and (val2, val1)
                mask1 = (abs(df[axis1] - val1) <= tolerance) & (abs(df[axis2] - val2) <= tolerance)
                mask2 = (abs(df[axis1] - val2) <= tolerance) & (abs(df[axis2] - val1) <= tolerance)

                combined_mask = mask1 | mask2
                matching_rows = df[combined_mask]

                for _, row in matching_rows.iterrows():
                    # Determine which orientation matches
                    if abs(row[axis1] - val1) <= tolerance and abs(row[axis2] - val2) <= tolerance:
                        target_axis1, target_axis2 = val1, val2
                    else:
                        target_axis1, target_axis2 = val2, val1

                    matches.append({
                        "refcode": row["refcode"],
                        "matched_axes": f"{axis1}_{axis2}",
                        f"{axis1}": row[axis1],
                        f"{axis2}": row[axis2],
                        f"target_{axis1}": target_axis1,
                        f"target_{axis2}": target_axis2,
                        f"diff_{axis1}": abs(row[axis1] - target_axis1),
                        f"diff_{axis2}": abs(row[axis2] - target_axis2),
                        "a": row["a"],
                        "b": row["b"],
                        "c": row["c"],
                        "SMILES": row["SMILES"],
                    })

        key = f"{val1}_{val2}"
        results[key] = pd.DataFrame(matches)

    return results


def save_matching_results(results: Dict[str, pd.DataFrame], output_dir: str = "output") -> None:
    """
    Save matching results to CSV files.

    Args:
        results: Dictionary mapping target values to DataFrames
        output_dir: Directory to save output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for target_values, df in results.items():
        filename = f"matched_{target_values}.csv"
        filepath = output_path / filename
        df.to_csv(filepath, index=False)


def main():
    """
    Main function to execute crystal lattice matching analysis.
    """
    # Load data
    data_path = Path("data/substituent/screening_refcode.parquet")
    df = pd.read_parquet(data_path)

    # Define target value pairs
    target_pairs = [(6.02, 7.23), (6.01, 13.2), (6.09, 9.74)]

    # Find matches for specific target values
    results = find_specific_value_matches(df, target_pairs, tolerance=0.25)

    # Display summary
    for target, matches_df in results.items():
        print(f"Target {target}: {len(matches_df)} matches found")

    # Save results
    save_matching_results(results)


if __name__ == "__main__":
    main()
