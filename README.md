## Replication Package for "Identifying Peer Effects in Networks with Unobserved 
Effort and Isolated Students

This folder contains the necessary files to replicate the results in the main 
paper and the online appendix. 

### R and Cpp functions
* Files `SourceR.R` and `SourceCpp.cpp` contain functions built in **R** and 
**cpp** to be used in other files.

### Monte Carlo Simulations
* `montecarlo.R` replicates our simulation study.

### Empirical Application
* This application uses data from Add Health, a program directed by Kathleen 
Mullan Harris and designed by J. Richard Udry, Peter S. Bearman, and Kathleen 
Mullan Harris at the University of North Carolina at Chapel Hill, and funded by 
Grant P01-HD31921 from the Eunice Kennedy Shriver National Institute of Child 
Health and Human Development, with cooperative funding from 23 other federal 
agencies and foundations. Information on how to obtain Add Health data files is 
available on the Add Health website (https://addhealth.cpc.unc.edu/).

* `0A_Inschool.do` extracts the part of the data set to be used from the Add 
Health data set.
* `0B_data.R` prepares the data set to be used.
* `1_exogenous_network.R` replicates the peer effect model estimation assuming 
that the network is exogenous
* `2_network_formation.R` estimates the network formation model.
* `3_endogenous_network.R` replicates the peer effect model estimation 
controlling for network endogeneity.
* `4_shocks.R` simulates shocks on alpha_s and c_s.

