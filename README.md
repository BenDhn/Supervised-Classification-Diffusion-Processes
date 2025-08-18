# Supervised Classification of Diffusion Processes

This repository presents a supervised classification framework for diffusion processes, combining rigorous mathematical modeling with simulations and machine learning techniques.

---

## Overview

This project includes:

- **Research Paper**  
  `Supervised_Classification_for_Diffusion_Processes.pdf` — A comprehensive document presenting the mathematical background, assumptions, and classification approach for stochastic differential equations (SDEs).

- **Model Implementations**  
  R scripts for simulating both homogeneous and non-homogeneous diffusion processes:
  - `MainSDEFunction_Homogeneous.R`
  - `UsefulFunctions_Homogeneous.R`
  - `MainSDEFunction_non_Homogeneous.R`
  - `UsefulFunctions_non_Homogeneous.R`
  - `SDEsimulation.R` — Example driver script for simulation

- **Predictions and Classification**  
  - `Predictions_non_homogeneous.Rmd` — R Markdown notebook for model application and prediction
  - `SDE_classification.ipynb` — Jupyter Notebook for interactive exploration and classification analysis

- **Other files**  
  - `.Rhistory`, `.RData`, `.DS_Store` — Optional environment-specific files

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/BenDhn/Supervised-Classification-Diffusion-Processes.git
cd Supervised-Classification-Diffusion-Processes
```

### 2.Set up your environment

R requirements: Install the following R packages if not already installed.

```r
install.packages(c("ggplot2", "dplyr", "caret", "deSolve", "knitr", "rmarkdown"))
```

Python requirements (for Jupyter Notebook):

```bash
pip install numpy pandas matplotlib scikit-learn jupyter
```

### 3. Run the simulations and classification

- To simulate and classify homogeneous or non-homogeneous SDEs, run the appropriate R script.
- To generate predictions and visuals, open and knit the Predictions_non_homogeneous.Rmd file.
- For a Python-based view, launch the SDE_classification.ipynb in Jupyter.
