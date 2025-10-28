# Extracted Matrices and Data

This directory contains extracted matrices and data from experimental studies for use in hiddenCOMET analysis.

## Directory Structure

```
extracted_matrices_and_data/
├── b_matrices/                    # B matrices from TGFB cases
│   ├── A549_TGFB_B_matrix.csv
│   ├── DU145_TGFB_B_matrix.csv
│   ├── MCF7_TGFB_B_matrix.csv
│   └── OVCA420_TGFB_B_matrix.csv
├── initial_emt_distributions/     # Initial EMT distributions from TGFB cases
│   ├── A549_TGFB_initial_distribution.csv
│   ├── DU145_TGFB_initial_distribution.csv
│   ├── MCF7_TGFB_initial_distribution.csv
│   └── OVCA420_TGFB_initial_distribution.csv
└── cell_cycle_data/              # Cell cycle data for EGF and TNF treatments
    ├── A549_EGF_cell_cycle.csv
    ├── A549_TNF_cell_cycle.csv
    ├── DU145_EGF_cell_cycle.csv
    ├── DU145_TNF_cell_cycle.csv
    ├── MCF7_EGF_cell_cycle.csv
    ├── MCF7_TNF_cell_cycle.csv
    ├── OVCA420_EGF_cell_cycle.csv
    └── OVCA420_TNF_cell_cycle.csv
```

## B Matrices (Emission Matrices)

### Description
B matrices represent the relationship between hidden EMT states and observed cell cycle phases. Each matrix shows P(cell_cycle_phase | EMT_state).

### Format
- **Rows**: EMT states (Epithelial, Hybrid, Mesenchymal)
- **Columns**: Cell cycle phases (G1_phase, S_phase, G2M_phase)
- **Values**: Probabilities (sum to 1 across each row)

### Usage
These matrices can be used directly in the hiddenCOMET app for:
- Main fitting workflow (Generator matrix estimation)
- Comparing EMT-cell cycle relationships across cell lines

## Initial EMT Distributions

### Description
Initial EMT state distributions at timepoint 0 for TGFB-treated cell lines.

### Format
- **Columns**: EMT states (Epithelial, Hybrid, Mesenchymal)
- **Values**: Probabilities (sum to 1)

### Usage
These distributions can be used as:
- Initial state distributions (p0) in the main fitting workflow
- Starting points for EMT trajectory modeling

## Cell Cycle Data

### Description
Time-series cell cycle phase data for EGF and TNF treatments across all cell lines.

### Format
- **Timepoints**: 0, 0.33, 1, 3, 7, 7.33, 8, 10 days
- **Columns**: G1_phase, S_phase, G2M_phase
- **Values**: Proportions (sum to 1 across each row)

### Usage
These datasets can be used for:
- Transition rate estimation
- Comparing cell cycle dynamics across treatments
- Input for observed process (Θ) in the main fitting workflow

## Data Sources

### TGFB Data
- **Source**: Real experimental data from TGFB-treated cell lines
- **Extraction**: B matrices and initial distributions derived from TGFB cell cycle data
- **Timepoints**: Full time series (0-10 days)

### EGF Data
- **Source**: EGF-treated cell lines
- **Purpose**: Cell cycle dynamics for transition rate estimation
- **Timepoints**: Full time series (0-10 days)

### TNF Data
- **Source**: TNF-treated cell lines
- **Purpose**: Cell cycle dynamics for transition rate estimation
- **Timepoints**: Full time series (0-10 days)

## Usage Examples

### For Main Fitting Workflow
1. Use B matrices from `b_matrices/` directory
2. Use initial distributions from `initial_emt_distributions/` directory
3. Use cell cycle data from `cell_cycle_data/` as observed process (Θ)

### For Transition Rate Estimation
1. Use EGF or TNF cell cycle data as observed process
2. Combine with appropriate B matrix and initial distribution
3. Estimate generator matrix (G) for treatment-specific dynamics

## Cell Line Characteristics

### A549 (Lung Cancer)
- **TGFB Response**: Strong EMT progression
- **Initial State**: Balanced epithelial-hybrid distribution
- **Cell Cycle**: Moderate G1 phase, high S phase

### DU145 (Prostate Cancer)
- **TGFB Response**: Moderate EMT with some reversal
- **Initial State**: Predominantly epithelial
- **Cell Cycle**: Balanced distribution

### MCF7 (Breast Cancer)
- **TGFB Response**: Starts mesenchymal, complex dynamics
- **Initial State**: Predominantly mesenchymal
- **Cell Cycle**: High G1 phase, low S phase

### OVCA420 (Ovarian Cancer)
- **TGFB Response**: Clear epithelial-to-mesenchymal transition
- **Initial State**: Predominantly epithelial
- **Cell Cycle**: Balanced distribution

## Contact

For questions about the analysis: annicenajafi27@gmail.com
