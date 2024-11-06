# Uncertainty-I2M2
Scripts for modeling uncertainty of the French WFD index I2M2

# Script 1
File: scripts/1_incertitudes_interoperator variability.R
This script takes the results from a preliminary inter-operator variability study and uses them to simulate 10,000 values for each site sampling event and each metric constitutive of the I2M2 index.
This script is also used to generate Table 1 and Figure S2.

# Script 2
File: scripts/2_incertitudes_index parameters.R
This script allows for the calculations of parameters used in the calculation of the I2M2, either as they were calculated during the development of the I2M2 (ie. Mondy et al. 2012), or calculated from randomized dataset (random sampling with replacement of the development dataset).

# Script 3
File: scripts/3_incertitudes_combinations of uncertainty sources.R
This script used the bootstrapped values of the parameters used in the calculation of the I2M2 (cf. script n°2), and/or Monte-Carlo simulations of the original metric values (cf. script n°1), before using these values in the bioindication process in order to simulate the uncertainty related to each of these uncertainty sources.
Note: this script takes at least 15 hours for all calculations. This time could definitely be shortened by doing some more matrix calculations.

# Script 4:
File: scripts/4_incertitudes_post analyses.R
This final script does some final data preparation before generating all remaining figures.
