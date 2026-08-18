import os
import numpy as np
import pandas as pd


# ===================== USER SETTINGS =====================

INPUT_FOLDER = r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated\xlsx"

OUTPUT_XLSX = (
    r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated"
    r"\FreqDist_NFE_AFR_AMR.xlsx"
)


# ===================== GROUP / POPULATION DEFINITIONS =====================

# Values used in "GroupMax FAF group"
GROUP_CODES = {
    "NFE": "nfe",
    "AFR": "afr",
    "AMR": "amr",
}

GROUP_COL = "GroupMax FAF group"


# Population-specific AF columns
FREQ_COLS = {
    "NFE": "Frequency European non-Finnish",
    "AFR": "Frequency African/African American",
    "AMR": "Frequency Admixed American",
}


# Population-specific allele count columns
AC_COLS = {
    "NFE": "Allele Count European (non-Finnish)",
    "AFR": "Allele Count African/African American",
    "AMR": "Allele Count Admixed American",
}


# ===================== ENRICHMENT COMPARISONS =====================
#
# For each designated ancestry:
#
# comparison passes if:
#
#     ratio >= 10
#
# OR
#
#     AF in the OTHER ancestry == 0
#
# BOTH ancestry comparisons must pass overall.

COMPARISONS = {

    "NFE": [
        (
            "European non-Finnish Frequency/African African American Frequency",
            "AFR"
        ),
        (
            "European non-Finnish Frequency/Admixed Frequency",
            "AMR"
        ),
    ],

    "AFR": [
        (
            "African African American Frequency/European non-Finnish Frequency",
            "NFE"
        ),
        (
            "African African American Frequency/Admixed Frequency",
            "AMR"
        ),
    ],

    "AMR": [
        (
            "Admixed Frequency/European non-Finnish Frequency",
            "NFE"
        ),
        (
            "Admixed Frequency/African African American Frequency",
            "AFR"
        ),
    ],
}


# ===================== CRITERIA =====================

POP_MIN = 0.00005
POP_MAX = 0.01

MIN_ALLELE_COUNT = 2

ENRICH_RATIO = 10


# ===================== AF BUCKETS =====================
#
# AF = 0 is counted separately.

BIN_EDGES = np.array([
    0.0,
    1e-5,
    1.25e-5,
    2.5e-5,
    3.75e-5,
    5e-5,
    6.25e-5,
    7.5e-5,
    8.75e-5,
    1e-4,
    1.25e-4,
    2.5e-4,
    3.75e-4,
    5e-4,
    6.25e-4,
    7.5e-4,
    8.75e-4,
    1e-3,
    1.25e-3,
    2.5e-3,
    3.75e-3,
    5e-3,
    6.25e-3,
    7.5e-3,
    8.75e-3,
    1e-2,
    np.inf
], dtype=float)


# ===================== HELPERS =====================

def make_bin_labels(edges):
    """
    Create readable AF bucket labels.

    AF = 0 is its own bucket.
    """

    labels = ["0"]

    for low, high in zip(edges[:-1], edges[1:]):

        if low == 0:
            labels.append(f">0-{high:g}")

        elif np.isinf(high):
            labels.append(f">={low:g}")

        else:
            labels.append(f"{low:g}-{high:g}")

    return labels


def convert_ratio(series):
    """
    Convert pre-calculated ratio column to numeric.

    #DIV/0! is allowed to become NaN because enrichment
    explicitly checks whether the denominator ancestry AF == 0.
    """

    series = series.replace({
        "#DIV/0!": np.nan,
        "#DIV/0": np.nan,
    })

    return pd.to_numeric(
        series,
        errors="coerce"
    )


def count_frequency_buckets(series):
    """
    Count AF distribution.

    AF = 0 is counted separately.
    Missing/non-numeric AF values are ignored.
    """

    values = pd.to_numeric(
        series,
        errors="coerce"
    ).dropna()

    zero_count = int(
        (values == 0).sum()
    )

    positive_values = values[
        values > 0
    ]

    hist_counts, _ = np.histogram(
        positive_values,
        bins=BIN_EDGES
    )

    counts = np.concatenate([
        [zero_count],
        hist_counts
    ])

    return counts


# ===================== ALL GROUPMAX =====================

def all_groupmax_mask(df, ancestry):
    """
    All variants assigned to the designated ancestry
    according to GroupMax FAF group.

    No AF or AC filtering.
    """

    group_value = GROUP_CODES[ancestry]

    return (
        df[GROUP_COL]
        .astype(str)
        .str.strip()
        .str.lower()
        == group_value
    )


# ===================== POPMAX =====================

def popmax_mask(df, ancestry):
    """
    PopMax analysis subset:

        GroupMax FAF group == designated ancestry
        AND
        0.00005 < designated ancestry AF < 0.01
        AND
        designated ancestry allele count >= 2
    """

    group_mask = all_groupmax_mask(
        df,
        ancestry
    )

    freq = pd.to_numeric(
        df[FREQ_COLS[ancestry]],
        errors="coerce"
    )

    allele_count = pd.to_numeric(
        df[AC_COLS[ancestry]],
        errors="coerce"
    )

    return (
        group_mask
        &
        (freq > POP_MIN)
        &
        (freq < POP_MAX)
        &
        (allele_count >= MIN_ALLELE_COUNT)
    )


# ===================== ENRICHED =====================

def enriched_mask(df, ancestry):
    """
    Enriched subset:

    First:
        must satisfy PopMax criteria.

    Then for EACH other ancestry:

        ratio >= 10

    OR

        AF in other ancestry == 0

    Both ancestry comparisons must pass overall.

    Example for NFE:

        GroupMax == nfe
        AND 0.00005 < NFE AF < 0.01
        AND NFE AC >= 2

        AND

        (
            NFE / AFR >= 10
            OR
            AFR AF == 0
        )

        AND

        (
            NFE / AMR >= 10
            OR
            AMR AF == 0
        )
    """

    mask = popmax_mask(
        df,
        ancestry
    )

    for ratio_col, other_ancestry in COMPARISONS[ancestry]:

        ratio = convert_ratio(
            df[ratio_col]
        )

        other_freq = pd.to_numeric(
            df[FREQ_COLS[other_ancestry]],
            errors="coerce"
        )

        comparison_pass = (
            (ratio >= ENRICH_RATIO)
            |
            (other_freq == 0)
        )

        mask = (
            mask
            &
            comparison_pass
        )

    return mask


# ===================== BUILD DISTRIBUTION =====================

def build_distribution(
    df,
    subset_type,
    labels
):
    """
    Output format:

        rows    = AF buckets
        columns = NFE / AFR / AMR
    """

    result = pd.DataFrame({
        "AF Bucket": labels
    })

    for ancestry in GROUP_CODES:

        # ---------- ALL ----------
        if subset_type == "All":

            mask = all_groupmax_mask(
                df,
                ancestry
            )

        # ---------- POPMAX ----------
        elif subset_type == "PopMax":

            mask = popmax_mask(
                df,
                ancestry
            )

        # ---------- ENRICHED ----------
        elif subset_type == "Enriched":

            mask = enriched_mask(
                df,
                ancestry
            )

        else:

            raise ValueError(
                f"Unknown subset type: {subset_type}"
            )


        subset = df[
            mask
        ]


        # Histogram uses the AF of the ancestry
        # that the variant belongs to
        counts = count_frequency_buckets(
            subset[
                FREQ_COLS[ancestry]
            ]
        )


        result[ancestry] = (
            counts.astype(int)
        )

    return result


# ===================== ADD TOTAL =====================

def add_total_row(table):
    """
    Add total count at bottom of each distribution table.
    """

    total_row = {
        "AF Bucket": "TOTAL"
    }

    for ancestry in GROUP_CODES:

        total_row[ancestry] = int(
            table[ancestry].sum()
        )

    return pd.concat(
        [
            table,
            pd.DataFrame([total_row])
        ],
        ignore_index=True
    )


# ===================== MAIN =====================

def build_frequency_distribution(
    input_folder,
    output_xlsx
):

    labels = make_bin_labels(
        BIN_EDGES
    )

    all_data = []

    print()
    print("Reading Excel files...")
    print()


    # ===================== READ FILES =====================

    for filename in sorted(
        os.listdir(input_folder)
    ):

        if not filename.lower().endswith(
            (".xlsx", ".xls")
        ):
            continue


        path = os.path.join(
            input_folder,
            filename
        )


        # Prevent output from being read again
        if os.path.abspath(path) == os.path.abspath(
            output_xlsx
        ):
            continue


        try:

            df = pd.read_excel(
                path
            )


            # ===================== REQUIRED COLUMNS =====================

            required_columns = [

                GROUP_COL,

                *FREQ_COLS.values(),

                *AC_COLS.values(),

                *[
                    ratio_col
                    for comparisons
                    in COMPARISONS.values()
                    for ratio_col, _
                    in comparisons
                ]
            ]


            # Remove duplicate names
            required_columns = list(
                dict.fromkeys(
                    required_columns
                )
            )


            missing_columns = [

                col
                for col in required_columns
                if col not in df.columns
            ]


            if missing_columns:

                print(
                    f"ERROR: {filename}"
                )

                print(
                    "Missing required columns:"
                )

                for col in missing_columns:

                    print(
                        f"    {col}"
                    )

                print()

                continue


            # Keep source file for traceability
            df["SourceFile"] = filename


            all_data.append(
                df
            )


            print(
                f"OK: {filename} "
                f"({len(df):,} rows)"
            )


        except Exception as e:

            print(
                f"ERROR reading "
                f"{filename}: {e}"
            )


    # ===================== CHECK =====================

    if not all_data:

        print()
        print(
            "No valid Excel files found."
        )

        return


    # ===================== COMBINE FILES =====================

    combined_df = pd.concat(
        all_data,
        ignore_index=True
    )


    print()
    print(
        f"Total rows loaded: "
        f"{len(combined_df):,}"
    )


    # ===================== DISTRIBUTIONS =====================

    all_distribution = (
        build_distribution(
            combined_df,
            "All",
            labels
        )
    )


    popmax_distribution = (
        build_distribution(
            combined_df,
            "PopMax",
            labels
        )
    )


    enriched_distribution = (
        build_distribution(
            combined_df,
            "Enriched",
            labels
        )
    )


    # ===================== ADD TOTAL ROWS =====================

    all_distribution = add_total_row(
        all_distribution
    )

    popmax_distribution = add_total_row(
        popmax_distribution
    )

    enriched_distribution = add_total_row(
        enriched_distribution
    )


    # ===================== PRINT SUMMARY =====================

    print()
    print(
        "========== SUMMARY =========="
    )


    for ancestry in GROUP_CODES:

        all_n = int(
            all_groupmax_mask(
                combined_df,
                ancestry
            ).sum()
        )

        popmax_n = int(
            popmax_mask(
                combined_df,
                ancestry
            ).sum()
        )

        enriched_n = int(
            enriched_mask(
                combined_df,
                ancestry
            ).sum()
        )


        print()
        print(ancestry)

        print(
            f"  GroupMax total: "
            f"{all_n:,}"
        )

        print(
            f"  PopMax subset:  "
            f"{popmax_n:,}"
        )

        print(
            f"  Enriched:       "
            f"{enriched_n:,}"
        )


    # ===================== SAVE EXCEL =====================

    output_folder = os.path.dirname(
        output_xlsx
    )


    if output_folder:

        os.makedirs(
            output_folder,
            exist_ok=True
        )


    with pd.ExcelWriter(
        output_xlsx,
        engine="openpyxl"
    ) as writer:


        all_distribution.to_excel(
            writer,
            sheet_name="All",
            index=False
        )


        popmax_distribution.to_excel(
            writer,
            sheet_name="PopMax",
            index=False
        )


        enriched_distribution.to_excel(
            writer,
            sheet_name="Enriched",
            index=False
        )


    print()
    print(
        f"Saved: {output_xlsx}"
    )


# ===================== RUN =====================

if __name__ == "__main__":

    build_frequency_distribution(
        INPUT_FOLDER,
        OUTPUT_XLSX
    )