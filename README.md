# hiddenCOMET

**Hidden Markov Models for EMT and Cellular Processes**

Access web application related to this project <a href="https://najafiannice.shinyapps.io/hiddenComet/">here</a>.


`hiddenCOMET` is an R package that provides functions for fitting continuous-time Markov chain (CTMC) models with hidden states, specifically designed for analyzing epithelial-mesenchymal transition (EMT) dynamics coupled with secondary cellular processes such as cell cycle progression.


The package Supports for any number of observed states with fixed 3-state EMT model

## Installation

```r
devtools::install_github("AnniceNajafi/hiddenCOMET")
```

Please follow the below steps to run the pipeline

### 1. Estimate Emission Matrix

```r
#Load your time-course data
Theta <- read.csv("observed_states_data.csv", row.names = 1)
P <- read.csv("hidden_states_trajectory.csv", row.names = 1)

#Estimate emission matrix using least squares
B_est <- estimate_B_ls(Theta, P)

#Or using maximum likelihood
B_est <- estimate_B_mle(Theta, P, 
                       observed_states = rownames(Theta),
                       hidden_states = c("E", "H", "M"))
```

### 2. Fit Generator Matrix

```r
#Define initial distribution
p0 <- c(0.8, 0.15, 0.05)
names(p0) <- c("E", "H", "M")

#Define timepoints
times <- c(0, 1, 2, 3, 5, 7)

#Fit generator matrix
result <- fit_G_given_p0(Theta = Theta, 
                        B = B_est, 
                        times = times, 
                        p0 = p0, 
                        states = c("E", "H", "M"))

#Access results
G_matrix <- result$G          #Generator matrix
P_trajectory <- result$P      #Hidden state trajectory
transition_rates <- result$par #Transition rate parameters
```

## Core Functions

### Model Fitting
- `fit_G_given_p0()`: Fit generator matrix using multi-start optimization
- `negloglik_G_p0fixed()`: Compute negative log-likelihood for CTMC model
- `build_G()`: Construct generator matrix with linear topology

### Emission Matrix Estimation
- `estimate_B_ls()`: Least squares estimation of emission matrix
- `estimate_B_mle()`: Maximum likelihood estimation of emission matrix

## Data Format

### Observed States Data (Theta)
```csv
timepoint,state1,state2,state3,state4
0,0.45,0.35,0.15,0.05
1,0.40,0.40,0.15,0.05
2,0.35,0.45,0.15,0.05
```

### Hidden States Trajectory (P)
```csv
timepoint,Epithelial,Hybrid,Mesenchymal
0,0.80,0.15,0.05
1,0.70,0.25,0.05
2,0.60,0.30,0.10
```

## Mathematical Background

The package implements a continuous-time hidden Markov model where:

- **Hidden states**: EMT states (Epithelial, Hybrid, Mesenchymal)
- **Observed states**: Any cellular process (e.g., cell cycle stages)
- **Generator matrix Q**: Defines transition rates between hidden states
- **Emission matrix B**: Maps hidden states to observed states probabilistically

The model assumes:
- Linear topology: E ↔ H ↔ M
- Continuous-time transitions with exponential waiting times
- Multinomial observation model

## Citation

If you use `hiddenCOMET` in your research, please cite:

```bibtex
@software{hiddenCOMET2025,
  title={hiddenCOMET: Hidden Markov Models for EMT and Cellular Processes},
  author={Annice Najafi},
  year={2024},
  url={https://github.com/AnniceNajafi/hiddenCOMET}
}
```

## License

MIT License - see LICENSE file for details.

## Contributing

This is an open source software. Contributions are welcome. Please feel free to submit a Pull Request if you have ideas to improve the pipeline.



