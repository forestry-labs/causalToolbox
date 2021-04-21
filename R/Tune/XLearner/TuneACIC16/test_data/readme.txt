This folder contains covariates, simulated treatment, and simulated response variables for the causal inference challenge in the 2016 Atlantic Causal Inference Conference. For each of 20 conditions, treatment and response data were simulated from real-world data corresponding to 4802 individuals and 58 covariates.

Files:
  x.csv - matrix of covariates; categorical variables are coded as A/B/C/..., binary variables as 0/1, and real numbers are left alone
  zy_##.csv - the twenty sets of treatment and response variables corresponding to various simulation settings; treatment is column "z" and response is column "y"


Submission:

Contest submissions should be in a csv file "est.csv" with the columns:
  est,ci_lower,ci_upper
where "est" is the treatment effect estimate for the effect of treatment on the treated, and "ci_lower"/"ci_upper" together give a 95% confidence/credible interval for your estimate. Rows in the file correspond to each of the 20 simulated data sets by number.

The file "runtime.txt" should contain a ball-park estimate of the computational time used in coming up with your estimates and a description of computer resources utilized.

Finally, if the method used can provide individual treatment effect estimates, these can be included in files "est_ind_XX.csv", where XX runs from 1 to 20. The format should correspond to "est.csv", but instead the rows corresponding to individuals.

Email a zip/tarball of your results to: vjd4@nyu.edu.

Thanks for participating!

ACIC2016 Organizing Comittee