# infectionAgeHIV
A lightweight version of the Multiple Biomarker Model for predicting HIV-1 infection times found in the R package biophybreak.

This package depends on rjags, which requires a separate installation of JAGS (see https://cran.r-project.org/web/packages/rjags/rjags.pdf for more info)

Workflow:
- format data into two dataframes, one for each patient's last negative HIV test (if known), first positive test, and ART start dates, and one for biomarker and sequence data
- use "find.infection.ages" to infer probability distributions for the amount of time between infection and diagnosis for each patient
