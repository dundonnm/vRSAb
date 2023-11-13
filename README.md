# vRSAb
SPM functions in MATLAB to perform variational representation similarity analysis using bootstrapped nulls

## Data structure
Place a contrast file (contrast.txt) and the vRSAb.m file in a parent directory  
Then place subject-specific data in different sub-directories  
![vRSAb1](https://github.com/dundonnm/vRSAb/assets/39175662/587bbfaa-585c-47a3-b67d-2e4fda29d5e6)

## Contrast file
A text file where columns describe the contrasts. E.g., the classic 2x2 ANOVA (2 main effects + interaction):  
1 1 -1  
-1 1 1  
1 -1 1  
-1 -1 -1  

## Inside data sub-directories
Place a text file for each run. Here there are six runs.  
Also place a design_good.txt file that specifies the design matrix concatenated to cover all runs
![vRSAb2](https://github.com/dundonnm/vRSAb/assets/39175662/425d6a5b-da6b-43ce-8e90-7528ae99b991)

## Citation
For the general logic of variational RSA please cite:  
Friston, K. J., Diedrichsen, J., Holmes, E., & Zeidman, P. (2019). Variational representational similarity analysis. _NeuroImage_, 201, 115986.

For our bootstrapped extension, please also cite:  
Marneweck, M., Gardner, C., Dundon, N. M., Smith, J., & Frey, S. H. (2023). Reorganization of sensorimotor representations of the intact limb after upper but not lower limb traumatic amputation. _NeuroImage: Clinical_, 39, 103499.

