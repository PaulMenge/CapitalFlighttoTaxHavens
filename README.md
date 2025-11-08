# Seminar Work on Capital Flight following the Arab Spring and Its Destinations

## Overview

This project investigates the dynamics of **capital flight around the Arab Spring** and the **destinations** that received those outflows. It is organized with a clear separation between **Code**, **Data**, and **Plots**.

**Initial goals (to refine together):**

* Describe how capital flight evolved in the years surrounding the Arab Spring.
* Identify primary destination countries, classified as financial centers and/or Tax Havens.
* Explore inequality and structure of capital flows following political and entailed oil shock.

## Repository Structure

```text
.
├── Code/
│   ├── notebooks/
│   │   ├── 01_data_overview.ipynb
│   │   ├── 02_cleaning_and_harmonization.ipynb
│   │   ├── 03_descriptive_trends.ipynb
│   │   ├── 04_destination_mapping.ipynb
│   │   └── 05_modeling_and_robustness.ipynb
│   ├── scripts/
│   │   ├── clean_data.py
│   │   ├── build_features.py
│   │   ├── analyze_capital_flight.py
│   │   └── export_plots.py
│   └── utils/
│       ├── io.py
│       ├── plotting.py
│       └── helpers.py
├── Data/
│   ├── raw/               # Unmodified datasets (not tracked if large/privileged; consider Git LFS)
│   ├── interim/           # Intermediate files (outputs of cleaning)
│   └── processed/         # Analysis-ready datasets
├── Plots/
│   ├── figures/           # Final figures used in the write-up
│   ├── exploratory/       # Quick charts from EDA
│   └── tables/            # LaTeX/CSV tables for paper/appendix
├── .gitignore
├── requirements.txt       # Python package requirements (to be finalized)
├── README.md              # You are here
└── LICENSE                # Chosen license (TBD)

## Data Catalog (to be completed)

Provide one row per dataset used.

| File | Location  | Source/URL | Coverage | Key Fields | Notes |
| ---- | --------- | ---------- | -------- | ---------- | ----- |
|      | Data/raw/ |            |          |            |       |

## Methods

* **Measurement of capital flight:** define metric(s) and construction logic.
* **Event window:** pre-/post- periods surrounding specific Arab Spring milestones (to be specified).
* **Destination identification:** mapping flows to receiving countries/centers.
* **Empirical approach:** descriptive trends, difference-in-differences/event study, robustness.
* **Sensitivity checks:** alternative metrics, alternative windows, country exclusions.
