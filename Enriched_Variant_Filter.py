import os
import pandas as pd


def filter_and_split_enriched_variants(
    input_folder: str,
    output_folder: str,
    group_column: str,
    clinvar_column: str,
    ancestry_config: dict,
    freq_min: float = 0.00005,
    freq_max: float = 0.01,
    ratio_threshold: float = 10,
    allele_count_min: int = 2,
):
    """
    For each Excel file:

    For each ancestry (NFE / AMR / AFR), identify ancestry-enriched variants.

    A variant is considered enriched in an ancestry if:

    1. GroupMax FAF group == that ancestry
    2. freq_min < ancestry AF < freq_max
    3. ancestry allele count > allele_count_min
    4. ancestry / comparison ancestry 1 AF ratio >= ratio_threshold
       OR comparison ancestry is absent (ratio is NaN)
    5. ancestry / comparison ancestry 2 AF ratio >= ratio_threshold
       OR comparison ancestry is absent (ratio is NaN)

    The variant must satisfy BOTH pairwise ancestry comparisons.

    Enriched variants are then split into:
        P/LP
        B/LB
        VUS
        Blank

    Outputs:
        3 ancestries × 4 ClinVar buckets = 12 Excel files
    """

    os.makedirs(output_folder, exist_ok=True)

    # -------------------------
    # ClinVar classification buckets
    # -------------------------

    P_LP = {
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
    }

    B_LB = {
        "Benign",
        "Likely benign",
        "Benign/Likely benign",
    }

    VUS = {
        "Uncertain significance",
    }

    # -------------------------
    # Storage
    # -------------------------

    buckets = {
        ancestry: {
            "P_LP": [],
            "B_LB": [],
            "VUS": [],
            "Blank": [],
        }
        for ancestry in ancestry_config
    }

    counts = {
        ancestry: {
            "enriched": 0,
            "P_LP": 0,
            "B_LB": 0,
            "VUS": 0,
            "Blank": 0,
        }
        for ancestry in ancestry_config
    }

    total_files = 0

    # -------------------------
    # Process Excel files
    # -------------------------

    for filename in os.listdir(input_folder):

        if not filename.lower().endswith((".xls", ".xlsx")):
            continue

        total_files += 1

        path = os.path.join(input_folder, filename)

        print(f"\nProcessing: {filename}")

        try:
            df = pd.read_excel(path)

            # Normalize GroupMax column
            group_series = (
                df[group_column]
                .fillna("")
                .astype(str)
                .str.lower()
                .str.strip()
            )

            # Normalize ClinVar column
            clin_series = (
                df[clinvar_column]
                .fillna("")
                .astype(str)
                .str.strip()
            )

            # -------------------------
            # Test each ancestry
            # -------------------------

            for label, config in ancestry_config.items():

                group_code = config["group_code"]
                freq_col = config["frequency"]
                ac_col = config["allele_count"]

                ratio_col_1 = config["ratio_1"]
                ratio_col_2 = config["ratio_2"]

                # Force numeric values
                freq = pd.to_numeric(
                    df[freq_col],
                    errors="coerce"
                )

                allele_count = pd.to_numeric(
                    df[ac_col],
                    errors="coerce"
                )

                ratio_1 = pd.to_numeric(
                    df[ratio_col_1],
                    errors="coerce"
                )

                ratio_2 = pd.to_numeric(
                    df[ratio_col_2],
                    errors="coerce"
                )

                # -------------------------
                # ENRICHMENT FILTER
                # -------------------------

                enriched_mask = (
                    # Must be PopMax ancestry
                    (group_series == group_code)

                    # Rare-variant AF window
                    & (freq > freq_min)
                    & (freq < freq_max)

                    # Recurrent variant
                    & (allele_count > allele_count_min)

                    # ≥10x vs ancestry #1 OR absent there
                    & (
                        (ratio_1 >= ratio_threshold)
                        | ratio_1.isna()
                    )

                    # ≥10x vs ancestry #2 OR absent there
                    & (
                        (ratio_2 >= ratio_threshold)
                        | ratio_2.isna()
                    )
                )

                df_sub = df.loc[enriched_mask].copy()

                if df_sub.empty:
                    print(f"  {label}: enriched=0")
                    continue

                # Add provenance columns
                df_sub.insert(
                    0,
                    "Enriched_Ancestry",
                    label
                )

                df_sub.insert(
                    0,
                    "Source_File",
                    filename
                )

                # ClinVar values corresponding to these rows
                clin_sub = clin_series.loc[df_sub.index]

                # -------------------------
                # ClinVar categories
                # -------------------------

                mask_p_lp = clin_sub.isin(P_LP)
                mask_b_lb = clin_sub.isin(B_LB)
                mask_vus = clin_sub.isin(VUS)
                mask_blank = clin_sub.eq("")

                buckets[label]["P_LP"].append(
                    df_sub.loc[mask_p_lp]
                )

                buckets[label]["B_LB"].append(
                    df_sub.loc[mask_b_lb]
                )

                buckets[label]["VUS"].append(
                    df_sub.loc[mask_vus]
                )

                buckets[label]["Blank"].append(
                    df_sub.loc[mask_blank]
                )

                # -------------------------
                # Counts
                # -------------------------

                counts[label]["enriched"] += len(df_sub)
                counts[label]["P_LP"] += int(mask_p_lp.sum())
                counts[label]["B_LB"] += int(mask_b_lb.sum())
                counts[label]["VUS"] += int(mask_vus.sum())
                counts[label]["Blank"] += int(mask_blank.sum())

                print(
                    f"  {label}: enriched={len(df_sub)} | "
                    f"P/LP={mask_p_lp.sum()} | "
                    f"B/LB={mask_b_lb.sum()} | "
                    f"VUS={mask_vus.sum()} | "
                    f"Blank={mask_blank.sum()}"
                )

        except Exception as e:
            print(f"  ERROR in {filename}: {e}")

    # -------------------------
    # Write outputs
    # -------------------------

    def write_bucket(
        ancestry_label: str,
        bucket_name: str,
        dfs: list[pd.DataFrame],
    ):

        out_path = os.path.join(
            output_folder,
            f"{ancestry_label}_{bucket_name}_enriched.xlsx"
        )

        if not dfs or all(d.empty for d in dfs):

            pd.DataFrame().to_excel(
                out_path,
                index=False
            )

            print(f"Saved (empty): {out_path}")
            return

        out_df = pd.concat(
            dfs,
            axis=0,
            ignore_index=True
        )

        out_df.to_excel(
            out_path,
            index=False
        )

        print(
            f"Saved: {out_path} "
            f"(rows={len(out_df)})"
        )

    for ancestry in ancestry_config:

        for bucket in [
            "P_LP",
            "B_LB",
            "VUS",
            "Blank",
        ]:

            write_bucket(
                ancestry,
                bucket,
                buckets[ancestry][bucket],
            )

    # -------------------------
    # Summary
    # -------------------------

    print("\n================ SUMMARY ================")

    print(f"Files processed: {total_files}")

    for ancestry in ancestry_config:

        c = counts[ancestry]

        print(
            f"{ancestry}: "
            f"enriched={c['enriched']} | "
            f"P/LP={c['P_LP']} | "
            f"B/LB={c['B_LB']} | "
            f"VUS={c['VUS']} | "
            f"Blank={c['Blank']}"
        )


# ============================================================
# YOUR SETTINGS
# ============================================================

file_path = "All Aggregated"

input_folder = rf"C:\Users\boxx_\Desktop\gnomAD\{file_path}\xlsx"

output_folder = (
    rf"C:\Users\boxx_\Desktop\gnomAD\{file_path}"
    rf"\outputs_enriched_by_ancestry"
)

group_column = "GroupMax FAF group"

clinvar_column = "ClinVar Germline Classification"


# ============================================================
# ANCESTRY-SPECIFIC COLUMNS
# ============================================================

ancestry_config = {

    "NFE": {

        "group_code": "nfe",

        "frequency":
            "Frequency European non-Finnish",

        "allele_count":
            "Allele Count European (non-Finnish)",

        "ratio_1":
            "European non-Finnish Frequency/Admixed Frequency",

        "ratio_2":
            "European non-Finnish Frequency/African African American Frequency",
    },


    "AMR": {

        "group_code": "amr",

        "frequency":
            "Frequency Admixed American",

        "allele_count":
            "Allele Count Admixed American",

        "ratio_1":
            "Admixed Frequency/European non-Finnish Frequency",

        "ratio_2":
            "Admixed Frequency/African African American Frequency",
    },


    "AFR": {

        "group_code": "afr",

        "frequency":
            "Frequency African/African American",

        "allele_count":
            "Allele Count African/African American",

        "ratio_1":
            "African African American Frequency/European non-Finnish Frequency",

        "ratio_2":
            "African African American Frequency/Admixed Frequency",
    },
}


# ============================================================
# RUN
# ============================================================

filter_and_split_enriched_variants(

    input_folder=input_folder,

    output_folder=output_folder,

    group_column=group_column,

    clinvar_column=clinvar_column,

    ancestry_config=ancestry_config,

    freq_min=0.00005,

    freq_max=0.01,

    ratio_threshold=10,

    allele_count_min=2,
)