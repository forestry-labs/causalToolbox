This folder contains covariates, simulated treatment, simulated response variables (including counterfactuals), and expected values for the causal inference challenge in the 2016 Atlantic Causal Inference Conference. For each of 20 conditions, treatment and response data were simulated from real-world data corresponding to 4802 individuals and 58 covariates.

Files:
  x.csv - matrix of covariates; categorical variables are coded as A/B/C/..., binary variables as 0/1, and real numbers are left alone
  zymu_##.csv - the twenty sets of treatment and response variables corresponding to various simulation settings; treatment is column "z", the observed response under the control is column "y0", the observed  response under treatment is "y1", the expected response under the control is "mu0", and the expected response under treatment is "mu1".

In this context, the ATT is calculated by:

  mean(mu1[z == 1] - mu0[z == 1]).

When differenced by an estimate, the result is collectively standardized by:

  sd(ifelse(z == 1, y1, y0))