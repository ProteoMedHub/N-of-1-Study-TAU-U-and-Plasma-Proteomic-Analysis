# N-of-1 HOMDD Studies TAU-U Analysis

## Overview
This repository contains comprehensive analysis code for studying protein expression patterns and treatment effects in patients undergoing TAU-U (Treatment as Usual - U) treatment. The project includes both proteomics data analysis and statistical evaluation of treatment outcomes.

## Project Structure
```
N-of-1-HOMDD-1-Studies-TAU-U-Analysis/
├── 0_PatA_stats.R        # Patient A statistical analysis
├── 0_PatK_stats.R        # Patient K statistical analysis
├── 1_PatA_protexpr.R     # Patient A protein expression analysis
├── 1_PatK_protexpr.R     # Patient K protein expression analysis
├── All_identified_proteins_PatA_filtered_norm.csv  # Patient A data
├── All_identified_proteins_PatK_filtered_norm.csv  # Patient K data
├── data.xlsx             # TAU-U treatment data
└── README.md            # This documentation
```

## Protein Expression Analysis

### Data Structure
The protein expression data is stored in CSV files with the following columns:
- `Accession`: Protein accession number
- `Peptide count`: Number of peptides identified
- `Unique peptides`: Number of unique peptides
- `Confidence score`: Score indicating protein identification confidence
- `Anova (p)`: ANOVA p-value for expression changes
- `q Value`: Adjusted p-value for multiple testing
- `Max fold change`: Maximum fold change in expression
- `Power`: Statistical power of the measurement
- `Highest mean condition`: Time point with highest mean expression
- `Lowest mean condition`: Time point with lowest mean expression
- `Mass`: Protein molecular mass
- `Description`: Protein description including gene name
- Multiple columns for expression measurements at different time points (W2, W4, W8, etc.)

### Analysis Scripts

#### 1. Statistical Analysis (0_*.R)
- Performs initial statistical analysis of protein expression data
- Conducts ANOVA and multiple testing corrections
- Identifies significant changes in protein expression

#### 2. Protein Expression Analysis (1_*.R)
- Analyzes protein expression patterns over time
- Creates visualizations including heatmaps
- Compares treatment vs placebo conditions

## TAU-U Treatment Effect Analysis

### Data Requirements
- Excel file containing treatment period data
- Columns: Treatment Period, PCS, MCS, BDI

### Statistical Methods
The analysis uses the `SingleCaseES` package to perform:
- Tau-U analysis (Parker's method)
- Linear regression for trend analysis
- Confidence interval calculations

### Visualizations
- Line plots for PCS, MCS, and BDI over time
- Faceted plots by treatment period
- Linear trend lines with confidence intervals

## Setup and Installation

### Required R Packages
```R
install.packages(c(
  "data.table",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "stringr",
  "readxl",
  "tidyverse",
  "scan",
  "kableExtra",
  "SingleCaseES",
  "pacman"
))
```

### Project Setup
1. Install all required R packages using the command above
2. Place the data files in their respective locations:
   - Protein expression data: `All_identified_proteins_PatA_filtered_norm.csv` and `All_identified_proteins_PatK_filtered_norm.csv`
   - TAU-U data: `data.xlsx`

## Usage

### Protein Expression Analysis
1. Ensure all required R packages are installed
2. Place the data files in the specified directory
3. Run the analysis scripts in sequence:
   ```R
   # First run the statistics scripts
   source("0_PatA_stats.R")
   source("0_PatK_stats.R")
   
   # Then run the protein expression analysis
   source("1_PatA_protexpr.R")
   source("1_PatK_protexpr.R")
   ```

### TAU-U Analysis
1. Ensure the Excel data file is properly formatted
2. Run the TAU-U analysis code:
   ```R
   # Load required packages
   pacman::p_load(tidyverse, scan, kableExtra, SingleCaseES)
   
   # Read data
   data <- read_excel("data.xlsx")
   
   # Run TAU-U analysis
   # (See code in README.md for complete analysis)
   ```

## Results Interpretation

### Protein Expression
- Significant changes in protein expression are identified using ANOVA and multiple testing corrections
- Heatmaps show expression patterns across different time points
- Treatment vs placebo comparisons highlight potential treatment effects

### TAU-U Outcomes
- Tau-U values indicate treatment effect size
- Confidence intervals provide statistical significance
- Visualizations show trends in PCS, MCS, and BDI scores

## Contact
For questions about the analysis or data, please contact the project maintainer.

