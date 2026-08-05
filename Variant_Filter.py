import os
import pandas as pd

def filter_and_split_all_ancestries(
    input_folder: str,
    output_folder: str,
    group_column: str,
    clinvar_column: str,
    freq_cols: dict,
    freq_min: float = 0.00005,
    freq_max: float = 0.01,
):
    """
    For each file:
      - For each ancestry (NFE/AMR/AFR):
          filter rows where GroupMax FAF group == ancestry_code
          AND freq_min < Frequency(ancestry) < freq_max
          then split into P/LP, B/LB, VUS, Blank and aggregate across files
    Writes 12 Excel files total: 3 ancestries * 4 ClinVar buckets.
    """

    os.makedirs(output_folder, exist_ok=True)

    # Buckets
    P_LP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
    B_LB = {"Benign", "Likely benign", "Benign/Likely benign"}
    VUS  = {"Uncertain significance"}

    # Your ancestry definitions
    ancestries = {
        "NFE": "nfe",  # non-Finnish European
        "AMR": "amr",  # Admixed American
        "AFR": "afr",  # African/African American
    }

    # Storage: buckets[ancestry_label][bucket_name] -> list of dfs
    buckets = {
        a: {"P_LP": [], "B_LB": [], "VUS": [], "Blank": []}
        for a in ancestries.keys()
    }

    total_files = 0
    counts = {a: {"matched": 0, "P_LP": 0, "B_LB": 0, "VUS": 0, "Blank": 0} for a in ancestries}

    for filename in os.listdir(input_folder):
        if not filename.lower().endswith((".xls", ".xlsx")):
            continue

        total_files += 1
        path = os.path.join(input_folder, filename)
        print(f"Processing: {filename}")

        try:
            df = pd.read_excel(path)

            # Normalize group col & clinvar col
            group_series = df[group_column].astype(str).str.lower().str.strip()
            clin_series  = df[clinvar_column].fillna("").astype(str).str.strip()

            for label, code in ancestries.items():
                freq_col = freq_cols[label]
                freq_num = pd.to_numeric(df[freq_col], errors="coerce")

                df_sub = df[
                    (group_series == code)
                    & (freq_num > freq_min)
                    & (freq_num < freq_max)
                ].copy()

                if df_sub.empty:
                    continue

                # provenance columns
                df_sub.insert(0, "Ancestry", label)
                df_sub.insert(0, "Source_File", filename)

                clin_sub = clin_series.loc[df_sub.index].fillna("").astype(str).str.strip()

                mask_p_lp  = clin_sub.isin(P_LP)
                mask_b_lb  = clin_sub.isin(B_LB)
                mask_vus   = clin_sub.isin(VUS)
                mask_blank = clin_sub.eq("")

                buckets[label]["P_LP"].append(df_sub.loc[mask_p_lp])
                buckets[label]["B_LB"].append(df_sub.loc[mask_b_lb])
                buckets[label]["VUS"].append(df_sub.loc[mask_vus])
                buckets[label]["Blank"].append(df_sub.loc[mask_blank])

                counts[label]["matched"] += len(df_sub)
                counts[label]["P_LP"] += int(mask_p_lp.sum())
                counts[label]["B_LB"] += int(mask_b_lb.sum())
                counts[label]["VUS"]  += int(mask_vus.sum())
                counts[label]["Blank"] += int(mask_blank.sum())

                print(
                    f"  {label}: matched={len(df_sub)} | "
                    f"P/LP={mask_p_lp.sum()} B/LB={mask_b_lb.sum()} "
                    f"VUS={mask_vus.sum()} Blank={mask_blank.sum()}"
                )

        except Exception as e:
            print(f"  ERROR in {filename}: {e}")

    # Write outputs (12 files)
    def write_bucket(ancestry_label: str, bucket_name: str, dfs: list[pd.DataFrame]):
        out_path = os.path.join(output_folder, f"{ancestry_label}_{bucket_name}_filtered.xlsx")
        if (not dfs) or all(d.empty for d in dfs):
            pd.DataFrame().to_excel(out_path, index=False)
            print(f"Saved (empty): {out_path}")
            return
        out_df = pd.concat(dfs, axis=0, ignore_index=True)
        out_df.to_excel(out_path, index=False)
        print(f"Saved: {out_path} (rows={len(out_df)})")

    for a in ancestries.keys():
        for b in ["P_LP", "B_LB", "VUS", "Blank"]:
            write_bucket(a, b, buckets[a][b])

    print("\nSummary")
    print(f"Files processed: {total_files}")
    for a in ancestries.keys():
        c = counts[a]
        print(
            f"{a}: matched={c['matched']} | "
            f"P/LP={c['P_LP']} B/LB={c['B_LB']} VUS={c['VUS']} Blank={c['Blank']}"
        )


# ---------------- YOUR SETTINGS ----------------
file_path = "All Aggregated"
input_folder = rf"C:\Users\boxx_\Desktop\gnomAD\{file_path}\xlsx"
output_folder = rf"C:\Users\boxx_\Desktop\gnomAD\{file_path}\outputs_split_by_ancestry"

group_column = "GroupMax FAF group"
clinvar_column = "ClinVar Germline Classification"

# Frequency columns per ancestry (edit if your headers differ)
freq_cols = {
    "NFE": "Frequency European non-Finnish",
    "AMR": "Frequency Admixed American",
    "AFR": "Frequency African/African American",
}

filter_and_split_all_ancestries(
    input_folder=input_folder,
    output_folder=output_folder,
    group_column=group_column,
    clinvar_column=clinvar_column,
    freq_cols=freq_cols,
    freq_min=0.00005,
    freq_max=0.01,
)
