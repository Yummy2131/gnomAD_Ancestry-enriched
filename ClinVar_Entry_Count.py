from pathlib import Path
import time
import re
import pandas as pd
from DrissionPage import ChromiumPage, ChromiumOptions

# ------------------- USER SETTINGS -------------------
INPUT_DIR   = Path(r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated\outputs_split_by_ancestry")   # <-- change this
OUTPUT_DIR  = Path(r"C:\Users\boxx_\Desktop\gnomAD\All Aggregated\outputs_split_by_ancestry\ClinVar") # <-- change this

ID_COL      = "gnomAD ID"         # Excel header must match
NEW_COL     = "ClinVar_Submissions"

GNOMAD_DATASET = "gnomad_r4"      # gnomAD v4.x
CHROME_PATH = r"/Applications/Google Chrome.app"
LOCAL_PORT  = 9202

SLEEP_BETWEEN_VARIANTS = 0.4
PAGE_TIMEOUT_SEC = 6.0

# Excel extensions
EXCEL_EXTS = {".xlsx", ".xlsm", ".xls"}  # .xls may or may not work depending on engine
# -----------------------------------------------------


def build_variant_url(gnomad_id: str) -> str:
    return f"https://gnomad.broadinstitute.org/variant/{gnomad_id}?dataset={GNOMAD_DATASET}"


def safe_get(page: ChromiumPage, url: str, timeout_sec: float) -> bool:
    try:
        page.get(url)
        start = time.time()
        while time.time() - start < timeout_sec:
            body = page.ele("tag:body")
            body_txt = (body.text or "").lower() if body else ""
            if "not found" in body_txt or "no results" in body_txt:
                return True
            if "clinvar" in body_txt or "submissions" in body_txt:
                return True
            time.sleep(0.2)
        return True
    except Exception:
        return False


def parse_clinvar_submissions(page: ChromiumPage) -> int:
    btn = page.ele('xpath://button[contains(., "See all") and contains(., "submissions")]', timeout=1)
    if btn:
        txt = (btn.text or "").strip()
        m = re.search(r"See all\s+(\d+)\s+submissions", txt)
        if m:
            return int(m.group(1))

    # fallback: look at a few candidate buttons only
    candidates = page.eles('xpath://button[contains(., "submissions")]', timeout=1) or []
    for b in candidates[:10]:
        txt = (b.text or "").strip()
        m = re.search(r"See all\s+(\d+)\s+submissions", txt)
        if m:
            return int(m.group(1))

    return 0


def get_clinvar_for_id(page: ChromiumPage, gnomad_id: str, cache: dict) -> int:
    gnomad_id = (gnomad_id or "").strip()
    if not gnomad_id or gnomad_id.lower() in ("nan", "none"):
        return 0
    if gnomad_id in cache:
        return cache[gnomad_id]

    url = build_variant_url(gnomad_id)
    ok = safe_get(page, url, PAGE_TIMEOUT_SEC)
    if not ok:
        cache[gnomad_id] = 0
        return 0

    n = parse_clinvar_submissions(page)
    cache[gnomad_id] = n
    time.sleep(SLEEP_BETWEEN_VARIANTS)
    return n


def process_workbook(page: ChromiumPage, xls_path: Path, out_path: Path, cache: dict) -> None:
    # Load all sheets as DataFrames
    sheets = pd.read_excel(xls_path, sheet_name=None)  # dict: {sheet_name: df}

    out_sheets = {}
    for sheet_name, df in sheets.items():
        if ID_COL not in df.columns:
            # leave sheet unchanged if no gnomAD ID col
            out_sheets[sheet_name] = df
            continue

        ids = df[ID_COL].astype(str).fillna("").tolist()
        clinvar_counts = [get_clinvar_for_id(page, vid, cache) for vid in ids]

        df2 = df.copy()
        # insert as first column; if it already exists, replace and move to front
        if NEW_COL in df2.columns:
            df2.drop(columns=[NEW_COL], inplace=True)
        df2.insert(0, NEW_COL, clinvar_counts)

        out_sheets[sheet_name] = df2

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df_out in out_sheets.items():
            df_out.to_excel(writer, sheet_name=str(sheet_name)[:31], index=False)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Start browser once
    co = ChromiumOptions().set_paths(browser_path=CHROME_PATH).set_local_port(LOCAL_PORT)
    page = ChromiumPage(co)

    # Cache across ALL files/sheets
    cache = {}

    files = [p for p in INPUT_DIR.rglob("*") if p.suffix.lower() in EXCEL_EXTS and not p.name.startswith("~$")]
    print(f"Found {len(files)} Excel files under: {INPUT_DIR}")

    for i, xls_path in enumerate(files, start=1):
        rel = xls_path.relative_to(INPUT_DIR)
        out_path = OUTPUT_DIR / rel  # preserve subfolders
        out_path = out_path.with_name(out_path.stem + "_with_clinvar" + out_path.suffix)

        print(f"[{i}/{len(files)}] Processing: {xls_path} -> {out_path}")
        try:
            process_workbook(page, xls_path, out_path, cache)
        except Exception as e:
            print(f"  ERROR: {xls_path} failed: {e}")

    print("Done.")


if __name__ == "__main__":
    main()
