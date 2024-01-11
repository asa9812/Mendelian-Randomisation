install.packages("remotes")                    # installing package
remotes::install_github("MRCIEU/TwoSampleMR")  # importing from github. 
library(TwoSampleMR)
library(ggplot2)
install.packages("MendelianRandomization")
library(MendelianRandomization)
# List available GWASs
ao <- available_outcomes() #
# Initiating a search to find the R library that includes the 'mr_egger' function.
search("R mr_egger function library")
help()
# Get instruments
#help("extract_instruments") #This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes. It then performs LD based clumping to return only independent significant associations.
#help("extract_instruments")
exposure_dat <- extract_instruments("ukb-b-19953") #THIS DATSET CONTAINS THE REGRESSION COEFFICIENTS.
head(exposure_dat)
summary(exposure_dat)
#help("clump_data")
clump_data(exposure_dat) # 19 variants of 458 variants removed due to LD with other variants or absence from LD reference panel. Note that extract_instruments performs clumping automatically.
# Get effects of instruments on outcome
#help("extract_outcome_data") #Supply the output from read_exposure_data() and all the SNPs therein will be queried against the requested outcomes in remote database using API.
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "ieu-b-39") #THIS DATSET CONTAINS THE REGRESSION COEFFICIENTS. ALL THE EXTRACTED DATA.help("extract_outcome_data")
#help("read_outcome_data") #Reads in outcome data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se.
#read_outcome_data(outcome_dat) # ? Note $SNP means accessing snps from the SNP column
head(outcome_dat)
summary(outcome_dat)
 
# Assuming SNP identifiers are in a column named 'SNP'
snps_list <- exposure_dat$SNP
head(snps_list)
 
 
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat) # Harmonise step.
#alternative harmonisaition from line 29 to line37.
ao <- available_outcomes()
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
 
# to harmonise do the following
dat <- harmonise_data(
    exposure_dat = bmi_exp_dat, 
    outcome_dat = chd_out_dat
# Perform MR
res <- mr(dat)
summary(res)
#extras
mr_singlesnp(dat)
mr_pleiotropy_test(dat) 
mr_heterogeneity(dat)
mr_steiger(dat)
res_single <- mr_singlesnp(dat)
# Leave one out analysis - It is possible to perform a leave-one-out analysis, where the MR is performed again but leaving out each SNP in turn, to identify if a single SNP is driving the association.
res_loo <- mr_leaveoneout(dat)
# plotting results
#Scatter plot - 
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
 
# A scatter plot is created for each exposure-outcome test, and stored in p1 as a list of plots. For example, to plot the first scatter plot:
p1[[1]]
 
# to see how many plots there are:
length(p1)
 
# Lines are drawn for each method used in mr(dat), the slope of the line corresponding to the estimated causal effect. To limit which lines are drawn, simply specify the desired methods, e.g. to only draw MR Egger and IVW:
res <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
p1 <- mr_scatter_plot(res, dat)
 
# It is possible to save this plot using the ggsave() function from the ggplot2 package, e.g. to save as a pdf
ggsave(p1[[1]], file = "filename.pdf", width = 7, height = 7)
 
# or save as PNG
ggsave(p1[[1]], file = "filename.png", width = 7, height = 7)


# Assume that 'exposure_dat' and 'outcome_dat' are the datasets containing summary statistics for BMI and diastolic blood pressure, respectively

# Harmonize the datasets
# harmonized_data <- harmonise_data(exposure_dat, outcome_dat)
# Load the MendelianRandomization package

# Assuming 'dat' is a data frame with columns for genetic variant (SNP) associations
# with the exposure (beta.exposure, se.exposure) and outcome (beta.outcome, se.outcome)
formatted_data <- mr_input(beta.exposure = dat$beta.exposure,
                           se.exposure = dat$se.exposure,
                           beta.outcome = dat$beta.outcome,
                           se.outcome = dat$se.outcome)

# Now apply the mr_egger function to the formatted data
mr_egger_result <- mr_egger(formatted_data)

# View the summary of the results
summary(mr_egger_result)


# IVW method with Q-test for heterogeneity
ivw_result <- ivw(dat, method="IVW")
summary(ivw_result)

# F-test for strength of instruments (example using a linear model, adjust as needed)
# Assuming 'instrument_strength' is a dataset or calculated value for the strength of the instrument
f_test_result <- summary(lm(instrument ~ exposure, data=instrument_strength))
f_test_result$fstatistic

# Over-identification test (Sargan test using MR-Base, example)
sargan_test_result <- mr_sargan(dat)
