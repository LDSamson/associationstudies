# associationstudies

This is a collection of functions that I used for longitudinal biomarker analyses in cohort studies. 

For more information on publications in which these functions are used, please visit [my personal website](https://ldsamson.github.io). 

The most important functions in this package are those that perform permutation testing between many different variables, using the permutation version of the spearman and wilcoxon tests as implemented in the coin package (see [here](https://doi.org/10.1198/000313006X118430) and [here](https://doi.org/10.18637/jss.v028.i08)) (coin::spearman_test and coin::wilcoxon_test, respectively). While permutation tests can be computational heavy, by using these non-parametric tests we do not have to make assumptions on the underlying distributions, hence making the tests more robust and reliable.

Other functions included in this package are functions facilitating prediction using random Forest, 
and functions that help to create publication-ready summary tables.

**Important: This is still a work in progress. Failures can be expected. More functions will be added in the future. Furthermore, I have to add more unit tests to test robustness of some functions**. 


To do: 

- Extend this readme file and description.
- Create more informative error message when data has <2 unique values in association study
- Create function that can test multiple pairs of associations, as implemented (here)[https://doi.org/10.1002/cti2.1374] (maybe merge with function association_study?)
- Merge the function association_study with association_study_long?
- Improve random forest function, implement unit tests
- Improve functions to create summary values for tables (geomean_conf)
- Create package vignettes
