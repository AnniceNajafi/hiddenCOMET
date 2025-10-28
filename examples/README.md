# hiddenCOMET Examples

This directory contains example files and usage demonstrations for the hiddenCOMET R package.

## Files

### Example Data Files
- **`example_B_matrix.csv`** - Example emission matrix B (observed states × hidden states)
- **`example_observed_process.csv`** - Example observed process data Θ (time series)
- **`example_p0_initial.csv`** - Example initial hidden state distribution p0
- **`sample_observed_states_data.csv`** - Sample observed states data for B matrix estimation

### Usage Examples
- **`example_usage.R`** - R script demonstrating how to use the hiddenCOMET functions
- **`app.R`** - Complete Shiny application for interactive analysis

## Quick Start

1. **Load the package:**
   ```r
   library(hiddenCOMET)
   ```

2. **Run the example script:**
   ```r
   source("examples/example_usage.R")
   ```

3. **Launch the Shiny app:**
   ```r
   shiny::runApp("examples/app.R")
   ```

## Data Format Requirements

### B Matrix (Emission Matrix)
- CSV format with observed states as rows and hidden states as columns
- Values should be probabilities (sum to 1 for each column)
- Example: 3 observed states (G1, S, G2M) × 3 hidden states (E, H, M)

### Observed Process Θ
- CSV format with timepoints as columns and observed states as rows
- Values should be fractions/proportions (sum to 1 for each column)
- First column should contain timepoint values

### Initial Distribution p0
- CSV format with hidden state names and initial probabilities
- Can be either two-column format (State, Value) or wide format
- Values should sum to 1

## Package Functions

The hiddenCOMET package provides the following main functions:

- `build_G()` - Build generator matrix with linear topology
- `fit_G_given_p0()` - Fit generator matrix given initial probabilities
- `estimate_B_ls()` - Least squares estimation of B matrix
- `estimate_B_mle()` - Maximum likelihood estimation of B matrix
- `proj_simplex()` - Project vectors onto probability simplex
- `check_matrix_conditioning()` - Check matrix numerical stability

For detailed documentation, see: `help(package = "hiddenCOMET")`
