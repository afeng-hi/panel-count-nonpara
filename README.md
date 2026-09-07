Simulation Studies
- **`Ser_Mean.R`** implements the proposed least-squares-based method for the simulation study in Section 4 of the main manuscript, where the recurrent-event count increments are generated from a two-point distribution.
- **`Ser_Mean_like.R`** implements the likelihood-based comparison method for the same simulation settings in Section 4 of the main manuscript.
- **`Ser_Mean_Poisson.R`** implements the proposed method under the nonhomogeneous Poisson process setting considered in Section 3.1 of the Supplementary Material.
- **`Ser_Mean_like_Poisson.R`** implements the corresponding likelihood-based method under the Poisson setting in Section 3.1 of the Supplementary Material.
- **`Ser_Mean_oracle.R`** implements the oracle estimator used as an additional benchmark in Section 3.1 of the Supplementary Material. 
- **`Ser_Mean_linear.R`** implements the linear method in Section 3.2 of the Supplementary Material.
- **`Ser_Mean_AFT.R`** implements the simulation study in Section 3.3 of the Supplementary Material, where the terminal event time is generated from a shifted accelerated failure time (AFT) model. 
- **`Summary.R`** summarizes the simulation output, including the calculation of bias and other reported performance measures, and produces the simulation figures.
- **`RE.R`** calculates the relative errors (REs) for the estimated nonparametric link function and baseline reversed mean function reported in the manuscript and Supplementary Material.
