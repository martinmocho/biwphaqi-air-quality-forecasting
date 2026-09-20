# Data

The datasets are **not included** in this repository. Download them from the original source and place them in this folder (they are git-ignored, so they won't be committed by accident). Please check the source's terms of use before redistributing any data.

| File expected by the script | Contents | Source |
|---|---|---|
| `climate_air_quality.csv` | Daily barangay-level concentrations of PM₂.₅, PM₁₀, NO₂, O₃, SO₂, and CO | [Project CCHAIN](https://doi.org/10.34740/KAGGLE/DS/4918229) |
| `worldpop_population.csv` | Annual barangay-level population (`pop_count_total`) | Project CCHAIN |

Pollutant sources within the dataset: NO₂, PM₁₀, and PM₂.₅ from CAMS; CO and SO₂ from MERRA-2; O₃ from ERA5. Project CCHAIN covers 12 Philippine cities across 2003 to 2022 and was produced by Thinking Machines Data Science, EpiMetrics, the Manila Observatory, and the Philippine Action for Community-led Shelter Initiatives, with support from the Lacuna Fund.

The script reads these files by name from the working directory, so run it from this folder (see the main README). It also writes an intermediate file, `monthly_adj_index_fixed.csv`, which is git-ignored.
