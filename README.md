# Filovirus Outbreak First Alerts: A Systematic Dataset (1976–2025)

## About This Repository

This repository contains the public data and analysis code accompanying the manuscript:

> **"First Alerts of Filovirus Outbreaks: A Systematic Review of 46 Ebola and Marburg Virus Disease Outbreaks, 1976–2025"**

We systematically reviewed grey and scientific literature to characterize the first alerts that triggered recognition of 46 filovirus (Ebola and Marburg virus) outbreaks reported between 1976 and 2025. The study examined 99 source documents across five databases to describe who raised early suspicion, what clinical and epidemiological signals were present, and how long detection and confirmation took.

### Key Findings

- The majority of first alerts (87%) came from physicians or hospital staff, not routine surveillance systems.
- No outbreak was detected by indicator-based surveillance alone.
- Time from first alert to laboratory confirmation has generally decreased over time, suggesting improved laboratory capacity.
- Timelines from the retrospectively identified first suspected case to outbreak recognition remained prolonged.
- Clinical suspicion by frontline healthcare workers and event-based surveillance continue to be the dominant — but often late — signals of filovirus outbreaks.

---

## Repository Contents

| File | Description |
|------|-------------|
| `Final_Extractions_All_public.csv` | Public extraction dataset: 99 source-document rows covering 46 unique outbreaks (CSV) |
| `Final_Extractions_All_public.xlsx` | Same dataset in Excel format |
| `references_clean.txt` | Cleaned bibliography of all 98 source documents included in the systematic review |
| `originalcode.r` | R script for all statistical analyses and figures in the manuscript |
| `figures/fig6_time_by_HCW.png` | Figure 6 — Time to declaration by healthcare worker involvement |
| `figures/fig7_time_by_signal.png` | Figure 7 — Time to declaration by first signal type |

---

## Dataset Description

### `Final_Extractions_All_public.csv`

Each row represents one source document. Multiple documents may cover the same outbreak (linked via `Outbreak_ID`). The dataset contains **99 rows** and **48 columns**.

Key fields include:

| Column | Description |
|--------|-------------|
| `Outbreak_ID` | Unique identifier for each outbreak (46 unique values) |
| `Virus_Reported` | Ebola or Marburg virus variant |
| `Outbreak_Start_Year` / `Outbreak_End_Year` | Year range of the outbreak |
| `Country`, `Region`, `Village` | Geographic location |
| `Source_Initial_Outbreak_Report` | Who first reported the outbreak |
| `Date_Initial_Outbreak_Report` | Date the initial report was made (YYYY-MM-DD where available) |
| `Key_Observations_Signals` | Narrative description of early clinical/epidemiological signals |
| `Challenges_Early_Detection_Reporting` | Barriers to early recognition documented in the source |
| `Total_Cases_Reported` | Final case count |
| `Total_Deaths_Reported` | Final death count |
| `CFR_Reported` | Case fatality ratio as reported |

Dates are formatted as `YYYY-MM-DD` where parseable; otherwise `NA`. Internal quality-assurance columns have been removed from this public release.

---

## R Analysis Code

`originalcode.r` reproduces all statistical analyses in the manuscript. It reads from `Extractions_R.csv` (the outbreak-level analytical dataset, available from the corresponding author on reasonable request) and produces all manuscript figures.

### Requirements

```r
install.packages(c("tidyverse", "ggplot2", "readr", "dplyr", "forcats", "stringr"))
```

### Usage

```r
# Place Extractions_R.csv in your Downloads folder, then:
source("originalcode.r")
```

---

## Citation

If you use this dataset or code, please cite:

> Feero M, et al. First Alerts of Filovirus Outbreaks: A Systematic Review of 46 Ebola and Marburg Virus Disease Outbreaks, 1976–2025. *[Journal]*. 2025.

---

## Contact

For questions about the dataset or analysis, please open a GitHub issue or contact the corresponding author.

---

## License

Data and code are released for academic and public health research use. Please cite the manuscript when using this material.
