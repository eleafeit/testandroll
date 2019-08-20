## Test & Roll: Profit Maximizing A/B Tests

# README File

#### Elea McDonnell Feit, eleafeit@gmail.com
#### Ron Berman, ronber@wharton.upenn.edu
#### Updated 20 August 2019

This folder contains the files necessary to replicate the results in Feit & Berman (2019) Test & Roll: Profit-Maximizing A/B Tests in R. A copy of the paper and the presentation are also included for reference. 

### Test & Roll Functions
`nn_functions.R` contains functions for computing and evaluating test & roll sample sizes for the normal-normal model and comparing them to null hypothesis tests and Thompson sampling. The functions in `nn_functions.R` will be included in the `testroll` R package with more complete documentation when it is released. 

### Website Example
`website.R` is a script for replicating the website example. First, it generates synthetic data and uses it to fit the model in `website_model.stan` using Stan via the `rstan` package. This illustrates how the meta-analysis was done; the data used in the website example in the paper is proprietary and can not be released. Second, it uses estimates from the model to find the optimal sample size for a test & roll using the functions in `nn_functions.R`. Note that the second part can be done without running the first part, because the parameter estimates reported in the paper have been hard-coded.  

`website_regret_sensitivity.R` produces the comparison between profit-maximizing test & roll experiments and Thompson sampling. 

### Display Advertising Example
`display.R` is a script for replicating the display advertising example.  It fits the model in `display_model.stan` using the data reported in Lewis and Rao (2015), which is in the file `display_LewisRao2015Retail.csv`. Then it uses the parameter estimates to design a profit-maximizing test & roll. 

### Catalog Example
`catalog.R` is similar to `website.R` and provides code for generating synthetic data similar to the catalog data, fitting the meta-analysis model (based on `catalog_model.stan`), and then computing optimal test & roll sample sizes. This example illustrates asymetric test design. 


