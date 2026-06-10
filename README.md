readme_content <- '# fastrenameR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fastrenameR)](https://CRAN.R-project.org/package=fastrenameR)
[![R-CMD-check](https://github.com/Serigne-Fallu/fastrenameR/workflows/R-CMD-check/badge.svg)](https://github.com/Serigne-Fallu/fastrenameR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**fastrenameR** is an R package for fast and flexible renaming of **Pathoplexus** FASTA sequence headers using metadata files. Designed for epidemiological and genomic workflows.

### New in v0.2.0: Date Standardization

| Input Format | Example | Output |
|-------------|---------|--------|
| Month-Day-Year | `5-15-2024` | `2024-05-15` |
| Month-Year | `5-2024` | `2024-05-15` |
| Decimal year | `2024.5` | `2024-07-02` |
| Year only | `2024` | `2024-07-01` |

## Features

- Fast processing of large FASTA files
- Pathoplexus compatible TSV metadata exports
- Automatic sanitisation of special characters
- Gzip support for compressed `.gz` FASTA files
- Flexible column selection by name or index
- Date standardization to ISO format for BEAST/BEAUti
- Date validation before and after conversion

## Installation

\`\`\`r
devtools::install_github("Serigne-Fallu/fastrenameR")
\`\`\`

## Quick Start

### 1. Explore your metadata

\`\`\`r
library(fastrenameR)
show_columns("mpox_metadata.tsv")
\`\`\`

### 2. Basic renaming

\`\`\`r
fastrenameR("mpox_metadata.tsv",
            "mpox_sequences.fasta",
            id_col = "accessionVersion",
            out_cols = c("geoLocCountry", "sampleCollectionDate"))
\`\`\`

**Before:** `>PQ123456`
**After:** `>USA_2024-05-15`

### 3. Rename with date standardization for BEAST

\`\`\`r
fastrenameR("metadata.tsv", "sequences.fasta",
            id_col = "accessionVersion",
            out_cols = c("geoLocCountry", "sampleCollectionDate"),
            standardize_dates = TRUE)
\`\`\`

### 4. Validate dates

\`\`\`r
dates <- c("2024-05-15", "5-15-2024", "2024.5", "2024", "invalid")
validate_dates(dates, report_format = "summary")
\`\`\`

## Real-world use cases

### MPOX sequences for BEAST analysis

\`\`\`r
library(fastrenameR)

show_columns("mpox_G1_metadata_2026-04-28.tsv")

fastrenameR("mpox_G1_metadata_2026-04-28.tsv",
            "mpox_G1_sequences.fasta",
            id_col = "accessionVersion",
            out_cols = c("geoLocCountry", "sampleCollectionDate", "lineage"),
            standardize_dates = TRUE)
# Output headers: >USA_2024-05-15_B.1.1.7
\`\`\`

### Compressed FASTA files

\`\`\`r
fastrenameR("metadata.tsv",
            "sequences.fasta.gz",
            id_col = "accessionVersion",
            out_cols = c("geoLocCountry", "sampleCollectionDate"))
\`\`\`

### Custom date parameters

\`\`\`r
fastrenameR("metadata.tsv", "sequences.fasta",
            id_col = "accessionVersion",
            out_cols = c("geoLocCountry", "sampleCollectionDate"),
            standardize_dates = TRUE,
            date_mid_month_day = 1,
            date_mid_year_day = 90)
\`\`\`

## Complete phylogenetic workflow

\`\`\`r
library(fastrenameR)

show_columns("project_metadata.tsv")

output_file <- fastrenameR("project_metadata.tsv",
                           "raw_sequences.fasta",
                           id_col = "accession",
                           out_cols = c("country", "collection_date", "strain"),
                           standardize_dates = TRUE)

validate_dates(names(output_file), report_format = "summary")
\`\`\`

## Troubleshooting

**Column not found**
Use `show_columns()` to see exact column names.

**Duplicate IDs**
Ensure `id_col` has unique values — last duplicate is kept.

**Missing sequences**
Check that all FASTA headers match IDs in your metadata.

## Performance

| File size | Time |
|-----------|------|
| 10,000 sequences | ~2 sec |
| 100,000 sequences | ~15 sec |
| 1,000,000 sequences | ~2 min |

## Dependencies

- **Required**: readr, stringr, lubridate
- **Optional**: Biostrings

## Citation

\`\`\`r
citation("fastrenameR")
\`\`\`

## License

MIT © Serigne Fallou Mbacke NGOM

## Support

- Issues: https://github.com/Serigne-Fallu/fastrenameR/issues
- Email: serignefalloumb.ngom@gmail.com
'

writeLines(readme_content, "README.md")
cat("README.md cree avec succes!\n")
cat("Lignes:", length(readLines("README.md")), "\n")