# Using MoTR to Probe Gender Agreement in Russian

This repository contains the code, data, and analysis for the project **Using MoTR to Probe Gender Agreement in Russian**, which involves statistical analyses on MoTR data and previously collected eye-tracking data. The project is organized into the following folders and files:

## Folder Structure

### Data and Resources

- **`data/`**  
  Contains the preprocessed MoTR reading measures used in the analysis.

- **`ref/`**  
  Includes references and supporting documents related to the project, including the eye-tracking data from   
  [Fuchs et al. (2025)](https://www.glossa-journal.org/article/id/11173/). You need to download the eye-tracking data and put it in this folder in order to run some of the analysis.

### Code and Models

- **`MoTR_ET_analysis.Rmd`** / **`MoTR_ET_analysis.nb.html`**  
  Main analysis script for the MoTR and eye-tracking data.  
  - `.Rmd`: Source file in R Markdown.  
  - `.nb.html`: Generated HTML report.

- **`model_comparison.Rmd`** / **`model_comparison.html`**  
  Script for Bayes Factor model comparison.  
  - `.Rmd`: Source file in R Markdown.  
  - `.html`: Generated HTML report.

- **`Power_analysis.Rmd`** / **`Power_analysis.html`**  
  Script and report for the power analysis.  
  - `.Rmd`: Source file in R Markdown.  
  - `.html`: Generated HTML report.

- **`stan/`**  
  Contains the Stan scripts used for Bayesian modeling in the main analysis.

### Reports and Outputs

- **`images/`**  
  Contains visualizations and images generated during the analysis.

- **`stats/`**  
  Stores statistics generated during the analysis. These include intermediate results for sanity checks or recycling in other parts of the analysis, such as in power analysis or model comparison.

