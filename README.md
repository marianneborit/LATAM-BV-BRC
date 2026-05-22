# Quantifying Latin America's Representation in the BV-BRC Genomic AMR Database

This repository contains the data and analysis code for an audit of how Latin
America is represented in the BV-BRC genomic antimicrobial resistance (AMR)
database, focusing on *Escherichia coli* ciprofloxacin resistance records.

The code reproduces every statistic in Table I of the paper, the regional
breakdowns, and all four figures, starting from the raw BV-BRC exports.

## Files

- `ecoli_fulldataset_-_BVBRC_genome_amr.csv.zip` — AMR phenotype table (E. coli x ciprofloxacin), zipped
- `BVBRC_genome.csv.zip` — genome metadata table (geography, year, host, etc.), zipped
- `amr_audit_reproduce.py` — reproduces Table I and the regional statistics
- `make_figures.py` — reproduces Figures 1-4
- `table1_geographic_distribution.csv`, `latam_country_counts.csv`,
  `top20_countries.csv`, `temporal_records_by_region.csv`,
  `missingness_by_region.csv` — derived tables produced by the audit script

## Requirements

Python 3.10 or newer, with pandas and matplotlib:

```
pip install pandas matplotlib
```

## Reproducing the results

First unzip the two data files so the scripts can read them:

```
unzip "ecoli_fulldataset_-_BVBRC_genome_amr.csv.zip"
unzip "BVBRC_genome.csv.zip"
```

Then run both scripts from the repository root.

Reproduce Table I and the headline numbers:

```
python amr_audit_reproduce.py \
    --amr "ecoli_fulldataset_-_BVBRC_genome_amr.csv" \
    --genome "BVBRC_genome.csv" \
    --csvout .
```

Reproduce the figures:

```
python make_figures.py \
    --amr "ecoli_fulldataset_-_BVBRC_genome_amr.csv" \
    --genome "BVBRC_genome.csv" \
    --outdir .
```

## Expected output

The audit script prints the following summary (Table I):

| Region          | Records | % of database | Countries | Most recent |
|-----------------|--------:|--------------:|----------:|------------:|
| Latin America   |      26 |         0.01% |         3 |        2004 |
| Rest of World   |   2,326 |         1.32% |        25 |        2019 |
| Country Unknown | 174,101 |        98.67% |         - |           - |
| Total           | 176,453 |          100% |        28 |        2019 |

The Latin American records come from Brazil (20), Mexico (4), and Uruguay (2).

The figure script writes four PNGs:

- `fig1_country_distribution.png` — top 20 contributing countries
- `fig2_resistance_rates.png` — ciprofloxacin resistance rate, LATAM vs rest of world
- `fig3_temporal.png` — annual and cumulative records over time, by region
- `fig4_missing_data.png` — metadata completeness heatmap by region

## How the two tables are joined

Each AMR phenotype record carries a `Genome ID`, but the geographic and clinical
metadata (country, collection year, host, isolation source) lives only in the
genome metadata table. The scripts left-join the phenotype table onto the genome
table on `Genome ID`, casting the key to a string first because the two exports
store it with different types. Records whose `Genome ID` has no matching metadata
row, or whose country field is blank, are counted as "Country Unknown."

## Data source

Both CSV files were downloaded from the BV-BRC database
(https://www.bv-brc.org/), filtered for *Escherichia coli* and ciprofloxacin.
BV-BRC is updated continuously, so re-downloading at a later date may give
slightly different record counts. To reproduce the exact numbers reported in the
paper, use the files included in this repository rather than a fresh download.

## Definition of Latin America

The 21-country region used throughout: Argentina, Bolivia, Brazil, Chile,
Colombia, Costa Rica, Cuba, Dominican Republic, Ecuador, El Salvador, Guatemala,
Haiti, Honduras, Mexico, Nicaragua, Panama, Paraguay, Peru, Puerto Rico,
Uruguay, and Venezuela.
