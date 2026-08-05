import os
import numpy as np
import pandas as pd


# ===================== USER SETTINGS =====================
INPUT_FOLDER  = r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated\xlsx"
OUTPUT_XLSX   = r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated\FreqDist_NFE_AMR.xlsx"

GROUP_COL     = "GroupMax FAF group"   # or "GroupMax FAF group"
GROUP_A       = "nfe"
GROUP_B       = "afr"

A_FREQ_COL    = "Frequency European non-Finnish"
B_FREQ_COL    = "Frequency African/African American"#"Frequency Admixed American"

# "ratio" columns used for ENRICHED definition:
A_OVER_B_RATIO_COL = "European non-Finnish Frequency/African African American Frequency" #"European non-Finnish Frequency/Admixed Frequency"
B_OVER_A_RATIO_COL = 'African African American Frequency/European non-Finnish Frequency' #"Admixed Frequency/European non-Finnish Frequency"

# criteria
POP_MIN = 0.0005
POP_MAX = 0.01
ENRICH_RATIO = 10

# Histogram bin edges (edit here if you want different resolution)
# Includes a 0 bin, then your old-style cut points up to 0.01
BIN_EDGES = [
    0.0,
    1e-5, 1.25e-5, 2.5e-5, 3.75e-5, 5e-5, 6.25e-5, 7.5e-5, 8.75e-5, 1e-4,
    1.25e-4, 2.5e-4, 3.75e-4, 5e-4,
    6.25e-4, 7.5e-4, 8.75e-4, 1e-3,
    1.25e-3, 2.5e-3, 3.75e-3, 5e-3,
    6.25e-3, 7.5e-3, 8.75e-3, 1e-2
]
# Make sure top edge includes POP_MAX (0.01)
if BIN_EDGES[-1] < POP_MAX:
    BIN_EDGES.append(POP_MAX)

BIN_EDGES = np.array(BIN_EDGES, dtype=float)


# ===================== HELPERS =====================
def _as_ratio(series: pd.Series) -> pd.Series:
    """
    Convert ratio column to numeric while preserving '#DIV/0!' as +inf.
    """
    s = series.copy()
    s = s.replace("#DIV/0!", np.inf)
    return pd.to_numeric(s, errors="coerce")

def _bin_labels(edges: np.ndarray) -> list[str]:
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == 0 and lo == 0.0:
            labels.append(f"0-{hi:g}")
        else:
            labels.append(f"{lo:g}-{hi:g}")
    return labels

def _hist_counts(values: pd.Series, edges: np.ndarray) -> np.ndarray:
    v = pd.to_numeric(values, errors="coerce")
    v = v.dropna().to_numpy()
    # clip to [0, 0.01) style range so outliers don't create surprises
    v = v[(v >= edges[0]) & (v < edges[-1])]
    counts, _ = np.histogram(v, bins=edges)
    return counts

def _mask_popmax(freq: pd.Series) -> pd.Series:
    f = pd.to_numeric(freq, errors="coerce")
    return (f > POP_MIN) & (f < POP_MAX)

def _mask_enriched(freq: pd.Series, ratio: pd.Series) -> pd.Series:
    pop = _mask_popmax(freq)
    r = _as_ratio(ratio)
    # enriched if ratio >= 10 OR ratio is +inf (DIV/0!)
    enr = (r >= ENRICH_RATIO) | np.isinf(r)
    return pop & enr


# ===================== MAIN =====================
def build_frequency_distribution(input_folder: str, output_xlsx: str) -> None:
    bins = BIN_EDGES
    labels = _bin_labels(bins)

    # row labels like your old output
    index_col = pd.DataFrame({"Bin": ["FileName", "Group", "Total #"] + labels})

    sheets = {
        "All": [],
        "popMax": [],
        "Enriched": [],
    }

    for fname in sorted(os.listdir(input_folder)):
        if not fname.lower().endswith((".xls", ".xlsx")):
            continue

        path = os.path.join(input_folder, fname)
        try:
            df = pd.read_excel(path)

            # ------- group A -------
            dfA = df[df[GROUP_COL] == GROUP_A].copy()
            a_freq = dfA[A_FREQ_COL]
            a_ratio = dfA[A_OVER_B_RATIO_COL] if A_OVER_B_RATIO_COL in dfA.columns else pd.Series([], dtype=object)

            # ------- group B -------
            dfB = df[df[GROUP_COL] == GROUP_B].copy()
            b_freq = dfB[B_FREQ_COL]
            b_ratio = dfB[B_OVER_A_RATIO_COL] if B_OVER_A_RATIO_COL in dfB.columns else pd.Series([], dtype=object)

            # ALL
            a_all_counts = _hist_counts(a_freq, bins)
            b_all_counts = _hist_counts(b_freq, bins)

            # popMax
            a_pop = dfA[_mask_popmax(a_freq)]
            b_pop = dfB[_mask_popmax(b_freq)]
            a_pop_counts = _hist_counts(a_pop[A_FREQ_COL], bins)
            b_pop_counts = _hist_counts(b_pop[B_FREQ_COL], bins)

            # ENRICHED (based on popMax)
            if A_OVER_B_RATIO_COL in dfA.columns:
                a_enr = dfA[_mask_enriched(a_freq, a_ratio)]
            else:
                a_enr = dfA.iloc[0:0]  # empty if ratio col missing
            if B_OVER_A_RATIO_COL in dfB.columns:
                b_enr = dfB[_mask_enriched(b_freq, b_ratio)]
            else:
                b_enr = dfB.iloc[0:0]

            a_enr_counts = _hist_counts(a_enr[A_FREQ_COL], bins) if len(a_enr) else np.zeros(len(labels), dtype=int)
            b_enr_counts = _hist_counts(b_enr[B_FREQ_COL], bins) if len(b_enr) else np.zeros(len(labels), dtype=int)

            # build output columns: for each file, two columns (A then B)
            def pack_col(group_name: str, total_n: int, counts: np.ndarray) -> list:
                return [fname, group_name, int(total_n)] + counts.astype(int).tolist()

            # totals are "how many values contributed to the histogram range"
            a_all_total = int(pd.to_numeric(a_freq, errors="coerce").notna().sum())
            b_all_total = int(pd.to_numeric(b_freq, errors="coerce").notna().sum())

            col_all_A = pd.DataFrame({f"{fname}__{GROUP_A}": pack_col(GROUP_A, a_all_total, a_all_counts)})
            col_all_B = pd.DataFrame({f"{fname}__{GROUP_B}": pack_col(GROUP_B, b_all_total, b_all_counts)})

            col_pop_A = pd.DataFrame({f"{fname}__{GROUP_A}": pack_col(GROUP_A, len(a_pop), a_pop_counts)})
            col_pop_B = pd.DataFrame({f"{fname}__{GROUP_B}": pack_col(GROUP_B, len(b_pop), b_pop_counts)})

            col_enr_A = pd.DataFrame({f"{fname}__{GROUP_A}": pack_col(GROUP_A, len(a_enr), a_enr_counts)})
            col_enr_B = pd.DataFrame({f"{fname}__{GROUP_B}": pack_col(GROUP_B, len(b_enr), b_enr_counts)})

            sheets["All"].extend([col_all_A, col_all_B])
            sheets["popMax"].extend([col_pop_A, col_pop_B])
            sheets["Enriched"].extend([col_enr_A, col_enr_B])

            print(f"OK: {fname}")

        except Exception as e:
            print(f"ERROR: {fname}: {e}")

    # write
    os.makedirs(os.path.dirname(output_xlsx), exist_ok=True)
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for sheet_name, cols in sheets.items():
            if cols:
                out = pd.concat([index_col] + cols, axis=1)
            else:
                out = index_col
            out.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Saved: {output_xlsx}")


if __name__ == "__main__":
    build_frequency_distribution(INPUT_FOLDER, OUTPUT_XLSX)
