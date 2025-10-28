# Quick Reference - Extracted Matrices and Data

## B Matrices (Emission Matrices)
**Purpose**: Show relationship between EMT states and cell cycle phases
**Format**: P(cell_cycle_phase | EMT_state)

| Cell Line | File | Key Characteristics |
|-----------|------|-------------------|
| A549 | `A549_TGFB_B_matrix.csv` | High S phase in epithelial, high G1 in mesenchymal |
| DU145 | `DU145_TGFB_B_matrix.csv` | Balanced distribution across states |
| MCF7 | `MCF7_TGFB_B_matrix.csv` | High G1 phase, low S phase |
| OVCA420 | `OVCA420_TGFB_B_matrix.csv` | High S phase in epithelial, high G1 in mesenchymal |

## Initial EMT Distributions
**Purpose**: Starting EMT state distribution at timepoint 0
**Format**: P(EMT_state | t=0)

| Cell Line | File | Initial State |
|-----------|------|--------------|
| A549 | `A549_TGFB_initial_distribution.csv` | Balanced epithelial-hybrid |
| DU145 | `DU145_TGFB_initial_distribution.csv` | Predominantly epithelial |
| MCF7 | `MCF7_TGFB_initial_distribution.csv` | Predominantly mesenchymal |
| OVCA420 | `OVCA420_TGFB_initial_distribution.csv` | Predominantly epithelial |

## Cell Cycle Data
**Purpose**: Time-series cell cycle data for transition rate estimation
**Format**: Proportions of G1, S, G2M phases over time

### EGF Treatment
- `A549_EGF_cell_cycle.csv`
- `DU145_EGF_cell_cycle.csv`
- `MCF7_EGF_cell_cycle.csv`
- `OVCA420_EGF_cell_cycle.csv`

### TNF Treatment
- `A549_TNF_cell_cycle.csv`
- `DU145_TNF_cell_cycle.csv`
- `MCF7_TNF_cell_cycle.csv`
- `OVCA420_TNF_cell_cycle.csv`

## Quick Usage

### For Main Fitting (Generator Matrix Estimation)
1. **B Matrix**: Use from `b_matrices/`
2. **Initial Distribution**: Use from `initial_emt_distributions/`
3. **Observed Process**: Use from `cell_cycle_data/`

### For Transition Rate Estimation
1. **Observed Process**: Use EGF or TNF data from `cell_cycle_data/`
2. **B Matrix**: Use appropriate TGFB B matrix
3. **Initial Distribution**: Use appropriate TGFB initial distribution

## Timepoints
All data includes: 0, 0.33, 1, 3, 7, 7.33, 8, 10 days
