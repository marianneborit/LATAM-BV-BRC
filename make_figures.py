#!/usr/bin/env python3
"""
Regenerate Figures 1-4 for the LACCI 2026 BV-BRC LATAM AMR audit, end-to-end,
directly from the two BV-BRC CSV exports. No intermediate files required.

Figure numbering matches the paper:
  fig1_country_distribution.png  -> Fig. 1  (Top 20 countries by AMR record count)
  fig2_resistance_rates.png      -> Fig. 2  (Ciprofloxacin resistance rate, LATAM vs RoW)
  fig3_temporal.png              -> Fig. 3  (Records over time + cumulative, by region)
  fig4_missing_data.png          -> Fig. 4  (Data-completeness heatmap by region)

INPUTS (same directory or via --amr/--genome):
  - ecoli_fulldataset_-_BVBRC_genome_amr.csv   (AMR phenotype table)
  - BVBRC_genome.csv                            (genome metadata table)

OUTPUT (--outdir, default ./figures)

This script is self-contained and deterministic: running it on the same inputs
reproduces the figures in the paper.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

LATAM = [
    "Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica",
    "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "Guatemala",
    "Haiti", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay",
    "Peru", "Puerto Rico", "Uruguay", "Venezuela",
]

LATAM_RED = "#D7263D"
ROW_BLUE = "#1f6fc4"


def load_merged(amr_path, genome_path):
    amr = pd.read_csv(amr_path, low_memory=False)
    genome = pd.read_csv(genome_path, low_memory=False)
    amr["Genome ID"] = amr["Genome ID"].astype(str).str.strip()
    genome["Genome ID"] = genome["Genome ID"].astype(str).str.strip()
    cols = ["Genome ID", "Isolation Country", "Collection Year",
            "Geographic Group", "Host Name", "Isolation Source"]
    g = genome[[c for c in cols if c in genome.columns]].drop_duplicates("Genome ID")
    merged = amr.merge(g, on="Genome ID", how="left")
    return merged


def region_series(country):
    def r(x):
        if pd.isna(x):
            return "Country Unknown"
        return "Latin America" if x in LATAM else "Rest of World"
    return country.map(r)


# Fig. 1 -- Top 20 countries by AMR record count
def fig1_country_distribution(merged, outdir):
    c = merged["Isolation Country"]
    counts = c.dropna().value_counts().head(20).sort_values(ascending=True)
    colors = [LATAM_RED if name in LATAM else ROW_BLUE for name in counts.index]

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.barh(counts.index, counts.values, color=colors)
    ax.set_title("Top 20 Countries by AMR Record Count\nE. coli Ciprofloxacin - BV-BRC Database",
                 fontsize=14, fontweight="bold")
    ax.set_xlabel("Number of AMR Records", fontsize=12)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=LATAM_RED, label="Latin America"),
                       Patch(facecolor=ROW_BLUE, label="Rest of World")],
              loc="lower right", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    p = os.path.join(outdir, "fig1_country_distribution.png")
    plt.savefig(p, dpi=150, bbox_inches="tight")
    plt.close()
    return p


# Fig. 2 -- Ciprofloxacin resistance rate, LATAM vs Rest of World
def fig2_resistance_rates(merged, outdir):
    region = region_series(merged["Isolation Country"])
    valid = merged["Resistant Phenotype"].isin(["Resistant", "Susceptible"])
    rates, ns = {}, {}
    for label in ["Latin America", "Rest of World"]:
        sub = merged[(region == label) & valid]
        ns[label] = len(sub)
        rates[label] = (sub["Resistant Phenotype"] == "Resistant").mean() * 100 if len(sub) else 0.0
    n_total_region = {lab: int((region == lab).sum()) for lab in ["Latin America", "Rest of World"]}

    fig, ax = plt.subplots(figsize=(8, 6))
    x = ["Latin America", "Rest of World"]
    y = [rates["Latin America"], rates["Rest of World"]]
    bars = ax.bar(x, y, color=[LATAM_RED, ROW_BLUE], width=0.55)
    ax.set_title("Ciprofloxacin Resistance Rate\nLATAM vs Rest of World Isolates",
                 fontsize=14, fontweight="bold")
    ax.set_ylabel("Ciprofloxacin Resistance Rate (%)", fontsize=12)
    ax.set_ylim(0, max(max(y) * 1.25, 1))
    for b, val in zip(bars, y):
        ax.text(b.get_x() + b.get_width() / 2, val + max(y) * 0.02 + 0.05,
                f"{val:.1f}%", ha="center", fontsize=13, fontweight="bold")
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels([f"Latin America\n(n={n_total_region['Latin America']})",
                        f"Rest of World\n(n={n_total_region['Rest of World']:,})"])
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.99, 0.97,
            f"Valid phenotype records: LATAM n={ns['Latin America']}, RoW n={ns['Rest of World']:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")
    plt.tight_layout()
    p = os.path.join(outdir, "fig2_resistance_rates.png")
    plt.savefig(p, dpi=150, bbox_inches="tight")
    plt.close()
    return p, ns, rates


# Fig. 3 -- Temporal coverage (annual + cumulative) by region
def fig3_temporal(merged, outdir):
    region = region_series(merged["Isolation Country"])
    df = merged.copy()
    df["__r"] = region
    df["__y"] = pd.to_numeric(df.get("Collection Year"), errors="coerce")
    df = df.dropna(subset=["__y"])
    df = df[(df["__y"] >= 1985) & (df["__y"] <= 2024)]

    years = np.arange(int(df["__y"].min()), int(df["__y"].max()) + 1)
    latam_by = df[df["__r"] == "Latin America"].groupby("__y").size().reindex(years, fill_value=0)
    row_by = df[df["__r"] == "Rest of World"].groupby("__y").size().reindex(years, fill_value=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Temporal Coverage of AMR Records by Region", fontsize=15, fontweight="bold")

    ax1.bar(years, row_by.values, color=ROW_BLUE, label="Rest of World")
    ax1.bar(years, latam_by.values, color=LATAM_RED, label="Latin America")
    ax1.set_title("Records Over Time\n(Rest of World vs Latin America)", fontsize=12)
    ax1.set_xlabel("Collection Year"); ax1.set_ylabel("Number of Records")
    ax1.legend(); ax1.grid(alpha=0.3)

    ax2.plot(years, np.cumsum(row_by.values), color=ROW_BLUE, lw=2.5, label="Rest of World")
    ax2.plot(years, np.cumsum(latam_by.values), color=LATAM_RED, lw=2.5, label="Latin America")
    latam_total = int(latam_by.sum())
    ax2.annotate(f"LATAM total:\n{latam_total} records",
                 xy=(years[-1], np.cumsum(latam_by.values)[-1]),
                 xytext=(years[-1] - 9, max(np.cumsum(row_by.values)) * 0.25),
                 fontsize=10, color=LATAM_RED,
                 arrowprops=dict(arrowstyle="->", color=LATAM_RED))
    ax2.set_title("Cumulative Records Over Time\nby Region", fontsize=12)
    ax2.set_xlabel("Year"); ax2.set_ylabel("Cumulative Records")
    ax2.legend(); ax2.grid(alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(outdir, "fig3_temporal.png")
    plt.savefig(p, dpi=150, bbox_inches="tight")
    plt.close()
    return p


# Fig. 4 -- Data-completeness heatmap by region
def fig4_missing_data(merged, outdir):
    region = region_series(merged["Isolation Country"])
    fields = ["Isolation Country", "Host Name", "Isolation Source",
              "Geographic Group", "Collection Year"]
    field_labels = ["Isolation\nCountry", "Host Name", "Isolation\nSource",
                    "Geographic\nGroup", "Collection\nYear"]
    rows = ["Latin America", "Rest of World"]
    mat = np.zeros((len(rows), len(fields)))
    for i, reg in enumerate(rows):
        sub = merged[region == reg]
        for j, f in enumerate(fields):
            mat[i, j] = (sub[f].isna().mean() * 100) if (f in sub and len(sub)) else 0.0

    cmap = LinearSegmentedColormap.from_list("gr", ["#1a7a3c", "#f7f7a8", "#b51330"])
    fig, ax = plt.subplots(figsize=(11, 4.2))
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(fields))); ax.set_xticklabels(field_labels, fontsize=10)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=11, fontweight="bold")
    ax.set_title("Data Completeness by Region\n(% of records with missing values per field)",
                 fontsize=13, fontweight="bold")
    for i in range(len(rows)):
        for j in range(len(fields)):
            v = mat[i, j]
            ax.text(j, i, f"{v:.0f}%", ha="center", va="center",
                    fontsize=12, fontweight="bold",
                    color="white" if (v > 70 or v < 8) else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("% Missing Data", fontsize=10)
    plt.tight_layout()
    p = os.path.join(outdir, "fig4_missing_data.png")
    plt.savefig(p, dpi=150, bbox_inches="tight")
    plt.close()
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amr", default="ecoli_fulldataset_-_BVBRC_genome_amr.csv")
    ap.add_argument("--genome", default="BVBRC_genome.csv")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    merged = load_merged(args.amr, args.genome)

    print("fig1:", fig1_country_distribution(merged, args.outdir))
    p2, ns, rates = fig2_resistance_rates(merged, args.outdir)
    print("fig2:", p2, "| valid-phenotype n:", ns,
          "| rates:", {k: round(v, 1) for k, v in rates.items()})
    print("fig3:", fig3_temporal(merged, args.outdir))
    print("fig4:", fig4_missing_data(merged, args.outdir))


if __name__ == "__main__":
    main()
