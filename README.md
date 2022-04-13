# mCP dye perturbation
This file contains the key function that developed in the 
Purified meta-Cresol Purple dye perturbation: How it influences spectrophotometric pH measurements Marine Chemistry [DOI:10.1016/j.marchem.2020.103849](https://www.sciencedirect.com/science/article/abs/pii/S0304420320301031?via%3Dihub)

The MATLAB code (dyeperturbation.m) included is the function to calculate the theoretical meta Cresol Purple (mCP) perturbation on measuring seawater and brackish water pH.The MATLAB function (dyeperturbation.m) can help with assessing the best S and pH for the dye stock and calculate the dye perturbation on water sample pH. This program is modified from CO2SYS (https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/).

# Instructions
Please see example_dyeperturbation.m

Input the estimated TA, DIC, S range and indictor properties. The function will output a matrix containing dye perturabion. Further analysis would the pH that deliver the samllest dye perturbation.

An detailed example of using this code can be found in the supplemental material of [DOI:10.1016/j.marchem.2020.103849](https://www.sciencedirect.com/science/article/abs/pii/S0304420320301031?via%3Dihub)

# Questions?
Please contact: xinyuli@udel.edu
