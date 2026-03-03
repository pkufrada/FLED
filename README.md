# EFIRM_fragment_decomposition
FLED_EFIRM_fragment_decomposition
# EFIRM NSCLC/KRAS MATLAB Analysis Scripts

This repository contains MATLAB scripts for EFIRM-based fragment decomposition and ROC analysis.
All scripts are written to load input files from the **current working directory**.

## Files
### MATLAB scripts
- `NSCLC_final_En.m` — NSCLC fragment decomposition (Short/Medium/Long) and outputs plots/CSVs
- `NSCLC_ROC_final_En.m` — ROC analysis for NSCLC
- `KRAS_data_emb_final_En.m` — KRAS analysis pipeline

### Input data (place in the same folder as the scripts)
- `20250917 All NSCLC Clinical .csv`
- `20250829 Plasma GIT-NEG.xlsx`
- `20250829 Saliva GIT-NEG.xlsx`

## Requirements
- MATLAB R2021a or later (recommended)
- Optimization Toolbox (required if `fmincon` is used)

## How to run
1. Download this repository as a ZIP (Code → Download ZIP), and unzip it.
2. Open MATLAB and set the **Current Folder** to the unzipped repository folder.
3. Run the scripts from the MATLAB Command Window:

```matlab
NSCLC_final_En
NSCLC_ROC_final_En
KRAS_data_emb_final_En
