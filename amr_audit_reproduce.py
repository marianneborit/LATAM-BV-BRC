#!/usr/bin/env python3
"""
Reproducible analysis for:
"Left Behind: Quantifying the Absence of Latin America in Global Genomic
 Antimicrobial Resistance Databases..."  (LACCI 2026)

This script regenerates the paper's headline statistics, Table I, and the
underlying data for Figures 1-5 directly from two BV-BRC CSV exports.

INPUTS (place in the same directory, or pass via --amr / --genome):
  1. ecoli_fulldataset_-_BVBRC_genome_amr.csv
       BV-BRC AMR phenotype table, E. coli x ciprofloxacin (~176k rows).
       Key columns: "Genome ID", "Antibiotic", "Resistant Phenotype".
  2. BVBRC_genome.csv
       BV-BRC genome metadata table for E. coli (~21k rows).
       Key columns: "Genome ID", "Isolation Country", "Collection Year",
                    "Geographic Group", "Host Name", "Isolation Source".

HOW THE TWO TABLES RELATE:
  Each AMR phenotype record carries a Genome ID. Geographic/clinical metadata
  lives only in the genome table. We left-join the phenotype table to the
  genome table on Genome ID (after standardizing the key to string) so that
  each phenotype record inherits its isolate's country, year, etc.

USAGE:
  python3 amr_audit_reproduce.py
  python3 amr_audit_reproduce.py --amr path/to/amr.csv --genome path/to/genome.csv --csvout out_dir

OUTPUT:
  - Prints all headline numbers and Table I to stdout.
  - Writes machine-readable CSVs (Table I, country counts, temporal series,
    missingness) to the chosen output directory for use as supplementary data.
"""

import argparse
import os
import sys
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# The 21-nation Latin America definition used in the paper (Methodology III-B).
LATAM = [
    "Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica",
    "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "Guatemala",
    "Haiti", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay",
    "Peru", "Puerto Rico", "Uruguay", "Venezuela",
]

# Metadata fields whose completeness is audited (Results IV-D / Fig. 5).
META_FIELDS = ["Isolation Country", "Collection Year", "Geographic Group",
               "Host Name", "Isolation Source"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def standardize_id(series):
    """Coerce Genome ID to a clean string key for joining.

    The two BV-BRC exports store Genome ID with different dtypes (one as float,
    one as string), so a naive merge fails. We cast both to stripped strings.
    """
    return series.astype(str).str.strip()


def pct(n, d):
    return (n / d * 100.0) if d else 0.0


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def load_and_join(amr_path, genome_path):
    amr = pd.read_csv(amr_path, low_memory=False)
    genome = pd.read_csv(genome_path, low_memory=False)

    amr["Genome ID"] = standardize_id(amr["Genome ID"])
    genome["Genome ID"] = standardize_id(genome["Genome ID"])

    # Defensive: confirm the phenotype table is the ciprofloxacin export.
    if "Antibiotic" in amr.columns:
        abx = amr["Antibiotic"].dropna().str.lower().unique()
        non_cip = [a for a in abx if a != "ciprofloxacin"]
        if non_cip:
            print(f"  [warn] phenotype table contains non-ciprofloxacin rows: {non_cip}")
            amr = amr[amr["Antibiotic"].str.lower() == "ciprofloxacin"].copy()

    # Keep one metadata row per genome, then left-join onto every phenotype row.
    meta_cols = ["Genome ID"] + [c for c in META_FIELDS if c in genome.columns]
    g = genome[meta_cols].drop_duplicates("Genome ID")
    merged = amr.merge(g, on="Genome ID", how="left")

    return amr, genome, merged


def region_of(country):
    if pd.isna(country):
        return "Country Unknown"
    return "Latin America" if country in LATAM else "Rest of World"


def run(amr, genome, merged, csvout=None):
    n_total = len(merged)
    country = merged["Isolation Country"] if "Isolation Country" in merged else pd.Series([pd.NA] * n_total)

    # --- Join transparency: why is country missing? -----------------------
    gid_set = set(genome["Genome ID"])
    matched = merged["Genome ID"].isin(gid_set)
    n_matched = int(matched.sum())
    n_unmatched = int((~matched).sum())
    has_country = country.notna()
    n_matched_blank = int((matched & ~has_country).sum())

    print("=" * 70)
    print("BV-BRC E. coli x Ciprofloxacin AMR audit -- reproduction")
    print("=" * 70)
    print(f"Total AMR phenotype records (cipro): {n_total:,}")
    print(f"Unique AMR Genome IDs:               {amr['Genome ID'].nunique():,}")
    print(f"Genome metadata rows available:      {genome['Genome ID'].nunique():,}")
    print()
    print("Why country is missing (decomposition of the ~98.6% 'unknown'):")
    print(f"  - No matching genome-metadata row:        {n_unmatched:,} ({pct(n_unmatched, n_total):.2f}%)")
    print(f"  - Matched but country field blank:        {n_matched_blank:,} ({pct(n_matched_blank, n_total):.2f}%)")
    print(f"  => Total without usable country:          {n_unmatched + n_matched_blank:,} "
          f"({pct(n_unmatched + n_matched_blank, n_total):.2f}%)")
    print(f"  Records WITH country information:         {int(has_country.sum()):,} "
          f"({pct(int(has_country.sum()), n_total):.2f}%)")

    # --- Region assignment -------------------------------------------------
    region = country.map(region_of)
    n_latam = int((region == "Latin America").sum())
    n_row = int((region == "Rest of World").sum())
    n_unknown = int((region == "Country Unknown").sum())

    latam_rows = merged[region == "Latin America"]
    row_rows = merged[region == "Rest of World"]

    def max_year(df):
        y = pd.to_numeric(df.get("Collection Year"), errors="coerce").dropna()
        return int(y.max()) if len(y) else None

    def n_countries(df):
        return int(df["Isolation Country"].nunique()) if "Isolation Country" in df else 0

    # --- Table I -----------------------------------------------------------
    table1 = pd.DataFrame([
        ["Latin America", n_latam, round(pct(n_latam, n_total), 2),
         n_countries(latam_rows), max_year(latam_rows)],
        ["Rest of World", n_row, round(pct(n_row, n_total), 2),
         n_countries(row_rows), max_year(row_rows)],
        ["Country Unknown", n_unknown, round(pct(n_unknown, n_total), 2),
         None, None],
        ["TOTAL", n_total, 100.0,
         n_countries(merged), max_year(merged)],
    ], columns=["Region", "Total Records", "% of Database",
                "Countries Represented", "Most Recent Record"])

    print("\n" + "-" * 70)
    print("TABLE I.  Geographic distribution of E. coli ciprofloxacin AMR records")
    print("-" * 70)
    print(table1.to_string(index=False))

    # --- LATAM detail (Results IV-B) --------------------------------------
    print("\nLATAM records by country:")
    latam_counts = latam_rows["Isolation Country"].value_counts()
    print(latam_counts.to_string())
    latam_years = pd.to_numeric(latam_rows.get("Collection Year"), errors="coerce").dropna()
    if len(latam_years):
        print(f"LATAM collection-year range: {int(latam_years.min())}-{int(latam_years.max())}")

    # --- Top-20 countries (Fig. 2) ----------------------------------------
    top20 = merged.loc[has_country, "Isolation Country"].value_counts().head(20)
    print("\nTop 20 contributing countries (Fig. 2 data):")
    print(top20.to_string())

    # --- Temporal series (Fig. 4) -----------------------------------------
    tmp = merged.copy()
    tmp["__region"] = region
    tmp["__year"] = pd.to_numeric(tmp.get("Collection Year"), errors="coerce")
    temporal = (tmp.dropna(subset=["__year"])
                   .groupby(["__year", "__region"]).size()
                   .unstack(fill_value=0).sort_index())

    # --- Missingness by region (Fig. 5) -----------------------------------
    miss_rows = []
    for label, df in [("Latin America", latam_rows), ("Rest of World", row_rows)]:
        rec = {"Region": label}
        for f in META_FIELDS:
            rec[f] = round(pct(df[f].isna().sum(), len(df)), 1) if (f in df and len(df)) else None
        miss_rows.append(rec)
    missingness = pd.DataFrame(miss_rows)
    print("\nMissing-metadata rate by region (% missing, Fig. 5 data):")
    print(missingness.to_string(index=False))

    # --- Write supplementary CSVs -----------------------------------------
    if csvout:
        os.makedirs(csvout, exist_ok=True)
        table1.to_csv(os.path.join(csvout, "table1_geographic_distribution.csv"), index=False)
        latam_counts.rename("records").to_csv(os.path.join(csvout, "latam_country_counts.csv"))
        top20.rename("records").to_csv(os.path.join(csvout, "top20_countries.csv"))
        temporal.to_csv(os.path.join(csvout, "temporal_records_by_region.csv"))
        missingness.to_csv(os.path.join(csvout, "missingness_by_region.csv"), index=False)
        print(f"\nSupplementary CSVs written to: {csvout}/")

    return table1


def main():
    p = argparse.ArgumentParser(description="Reproduce the BV-BRC LATAM AMR audit.")
    p.add_argument("--amr", default="ecoli_fulldataset_-_BVBRC_genome_amr.csv",
                   help="Path to the AMR phenotype CSV.")
    p.add_argument("--genome", default="BVBRC_genome.csv",
                   help="Path to the genome metadata CSV.")
    p.add_argument("--csvout", default=None,
                   help="Optional directory to write supplementary CSVs.")
    args = p.parse_args()

    for path in (args.amr, args.genome):
        if not os.path.exists(path):
            sys.exit(f"ERROR: input file not found: {path}")

    amr, genome, merged = load_and_join(args.amr, args.genome)
    run(amr, genome, merged, csvout=args.csvout)


if __name__ == "__main__":
    main()
