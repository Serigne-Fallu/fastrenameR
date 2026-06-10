# fastrenameR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fastrenameR)](https://CRAN.R-project.org/package=fastrenameR)
[![R-CMD-check](https://github.com/Serigne-Fallu/fastrenameR/workflows/R-CMD-check/badge.svg)](https://github.com/Serigne-Fallu/fastrenameR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**fastrenameR** is an R package for fast and flexible renaming of **Pathoplexus** FASTA sequence headers using metadata files. It was designed specifically for epidemiological and genomic workflows where sequence headers need to be standardised with metadata fields like:

- 🌍 Country of collection
- 📅 Collection date
- 🧬 Lineage or clade
- 🔢 Accession numbers
- 📊 Any other metadata column

### New in v0.2.0: Date Standardization

Convert various date formats to ISO standard (YYYY-MM-DD) for BEAST/BEAUti compatibility:

| Input Format | Example | Output |
|-------------|---------|--------|
| Month-Day-Year | `5-15-2024` | `2024-05-15` |
| Month-Year | `5-2024` | `2024-05-15` |
| Decimal year | `2024.5` | `2024-07-02` |
| Year only | `2024` | `2024-07-01` |

## Features

- 🚀 **Fast processing** - Handles large FASTA files efficiently
- 📊 **Pathoplexus compatible** - Works directly with Pathoplexus TSV metadata exports
- 🔄 **Automatic sanitisation** - Converts spaces and punctuation to safe hyphens
- 📦 **Gzip support** - Reads compressed `.gz` FASTA files directly
- 🎯 **Flexible column selection** - Use column names or indices
- 📝 **Informative output** - See exactly which sequences were renamed
- 📅 **Date standardization** - Convert dates to ISO format for phylogenetic software
- ✅ **Date validation** - Check date formats before and after conversion

## Installation

### From GitHub (recommended)

```r
# Using devtools
devtools::install_github("Serigne-Fallu/fastrenameR")

# Using pak (faster)
pak::pak("Serigne-Fallu/fastrenameR")

# Using remotes
remotes::install_github("Serigne-Fallu/fastrenameR")