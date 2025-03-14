# associationstudies

This is a collection of functions that I used for longitudinal biomarker analyses in cohort studies. 

For more information on publications in which these functions are used, please visit [my personal website](https://ldsamson.github.io). 

The most important functions in this package are those that perform permutation testing between many different variables, using the permutation version of the spearman and wilcoxon tests as implemented in the coin package (see [here](https://doi.org/10.1198/000313006X118430) and [here](https://doi.org/10.18637/jss.v028.i08)) (coin::spearman_test and coin::wilcoxon_test, respectively). While permutation tests can be computational heavy, by using these non-parametric tests we do not have to make assumptions on the underlying distributions, hence making the tests more robust and reliable.

Other functions included in this package are functions facilitating prediction analyses using the random Forest algorithm, 
and functions that help to create publication-ready summary tables.
