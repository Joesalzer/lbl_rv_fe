# Run line-by-line fixed effects models
Usage: Rscript lblfe_fit.R \<inputCSVPath\> \<outputDir\>

Where <inputCSVPath> is the path to a csv file. The first four columns must be named "lineID", "timeID", "date", and "rv_dirty". Any additional columns are treated as covariates, and each line will be assigned a unique slope for each covariate.

Requires rv_helper.R in the working directory.

This script fits both OLS and IRLS versions of the model. Results are saved in the <outputDir> with the base name of the csv file, as <fileName>__fit=OLS.rds and <fileName>__fit=IRLS.rds.

# Run line-by-line fixed effects models, with CV
Usage: Rscript lblfe_cv_fit.R \<inputCSVPath\> 
<outputDir\>

This script fits both OLS and IRLS versions of the model. Results are saved in the <outputDir> with the base name of the csv file, as <fileName>__fit=OLS.rds and <fileName>__fit=OLS__predictions.csv.
