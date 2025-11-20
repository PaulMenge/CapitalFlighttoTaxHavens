# Seminar Work on Capital Flight following the Arab Spring

This project investigates the dynamics of **capital flight around the Arab Spring** and the **destinations** that received those outflows. It is organized with a clear separation between **Code**, **Data**, and **Plots**.

* Describe how capital flight evolved in the years surrounding the Arab Spring.
* Identify primary destination countries, classified as financial centers and/or Tax Havens.
* Explore inequality and structure of capital flows following political and entailed oil shock.

## Repository Structure

```text
├── Code/
│ ├── BigEventStudycountry-block bootstrap CIs.R
│ ├── DDD_with_16+-qartals.R
│ ├── DDDnew.R
│ ├── EventStudyAll.R
│ ├── EventStudyforTH10.R
│ ├── newnewDDD.R
│ ├── newnewDDDnonEU.R
│ ├── non-EU-TH10.R
│ └── non-EU.R
├── Plots/
│ ├── 10QES.png
│ ├── 10QESNonEU.png
│ ├── Contributions.png
│ ├── ContributionsNonEU.png
│ ├── Did2.png
│ ├── DidNonEU.png
│ ├── JackknifeNonEU.png
│ ├── JackknifeRobust.png
│ ├── PlaceboNonEU.png
│ ├── PlaceboShocks.png
│ ├── PreTrend.png
│ ├── ReContributions.png
│ ├── ReContributionsNonEU.png
│ ├── Window.png
│ └── WindowNonEU.png
├── Data/
│   └── BIS Locational Banking Statistic Links.md
├── README.md
└── LICENSE                #(TBD)
```

### Folder Purpose & Conventions

* **Code/**: R scripts implementing event studies (DID/DDD), robustness checks, and figure generation/export.

  * Current scripts:

    * `BigEventStudycountry-block bootstrap CIs.R`
    * `DDD_with_16+-qartals.R`
    * `DDDnew.R`
    * `EventStudyAll.R`
    * `EventStudyforTH10.R`
    * `newnewDDD.R`
    * `newnewDDDnonEU.R`
    * `non-EU-TH10.R`
    * `non-EU.R`
* **Plots/**: PNG outputs produced by the analysis (event-study figures, placebo tests, jackknife robustness, contributions, and window visualizations).

*Data Catalog
Should be downloaded following the links in the Data folder and placed there.
| File                          | Location  | Source/URL                                                               | Coverage                                                                      | Key Fields                                                                                                    | Notes                                                                                                             |
| ----------------------------- | --------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| BIS Locational Banking Statistic Links | Data/ | Bank for International Settlements (BIS) – Locational Banking Statistics | Release as of Jul 2025; quarterly cross-border banking positions by residence | reporting_country, counterparty_country, quarter, position_type (assets/liabilities), currency (USD), value | **Single source for all code.** Downloaded/ingested Jul 2025;|

---

## Methods

* **Measurement of capital flight:** define metric(s) and construction logic.
* **Event window:** pre-/post- periods surrounding specific Arab Spring timeline.
* **Destination identification:** mapping flows to receiving countries/centers.
* **Empirical approach:** descriptive trends, difference-in-differences/event study, robustness.
* **Sensitivity checks:** alternative metrics, alternative windows, country exclusions.
