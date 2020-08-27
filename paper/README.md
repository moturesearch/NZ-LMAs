# Delineating functional labour market areas with estimable classification stabilities

This folder contains the source files for the paper describing our LMA delineation method.

## Workflow

We generate data, figures and tables by running the R scripts `code/data.R` and `code/analysis.R` (in that order).
We run each script in a fresh `paper.Rproj` instance within [RStudio](https://www.rstudio.com/).
`logs/` contains the session information from when we ran each script.

### Data

The table below describes the files in `data/`.
Confidentiality rules prevent us from sharing the raw data from which we derive these files.

#### Inputs (raw data)

`TTWA_flows_<yy>.dta` (excluded from this repository) contains inter-area unit commuter counts in Census year `<yy>`.
It contains three variables:

* `ur_au13`: Area unit containing usual residence
* `wp_au13`: Area unit containing workplace
* `flow`: Number of commuters with usual residence in `ur_au13` and workplace in `wp_au13`

#### Outputs

File | Description
--- | ---
`allocations.csv` | Crosswalk between area units (`au13`), LMAs (`lma`), and sub-LMAs (`sub_lma`) in each census year (`year`). Also includes LMA and sub-LMA classification stabilities.
`distance-thresholds.csv` | Sum of commuting flows included in our cleaned data (our "sample"), aggregated by estimated commute distance (in 50km bins) and census year.
`ensembles-main.csv` | Communities to which area units are assigned in each of 250 runs of the Louvain algorithm, by census year.
`ensembles-sub.csv` | Sub-communities to which area units are assigned in each of 250 runs of the Louvain algorithm, by census year.
`geographies.csv` | Crosswalk between area units, regional councils (`rc13`) and territorial authorities (`ta13`), derived from the 2013 Annual Areas file available [here](http://archive.stats.govt.nz/browse_for_stats/Maps_and_geography/Geographic-areas/geographic-area-files.aspx) (retrieved July 8, 2019).
`inter-sub-lma-flows.csv` | Table of inter-sub-LMA flows by origin sub-LMA, destination sub-LMA, and census year. Excludes flows with raw values below six.
`national-coverage.csv` | Table of total flows (`flow_total`), sample flows (`flow_sample`), flows with residential or workplace addresses belonging to unknown, outside-TA, employee-less or landless area units (`flow_bad_au13`), and flows with commute distances greater than 150km (`flow_long_commute`), aggregated by census year.
`populations.csv` | Table of sample employed resident (`employed_resident`), sample employee (`employee`), sample resident employee (`resident_employee`), and total employed resident (`employed_resident_total`) populations by area unit and census year.
`regional-coverage.csv` | Table of total and sample outflows, aggregated by regional council and census year.
`residual-sub-lma-flows.csv` | Table of inflows (`residual_inflow`) and outflows (`residual_outflow`) excluded from `inter-sub-lma-flows.csv`, aggregated by sub-LMA and census year.

### Dependencies

We use several R packages, identified at the start of the scripts in `code/` and in the log files in `logs/`.
All required packages are available on [CRAN](https://cran.r-project.org/).
