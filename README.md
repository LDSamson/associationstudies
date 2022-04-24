# associationstudies

Custom functions that are useful to perform association studies and prediction analyses in biomedical data, in which associations between many variables need to be tested.

I used these functions during my PhD. For more information on publications 
in which these functions are used, please visit [my personal website](https://ldsamson.github.io). 

This package relies heavily on the use of the permutation version of the spearman and wilcoxon test from the coin package (see [here](https://doi.org/10.1198/000313006X118430) and [here](https://doi.org/10.18637/jss.v028.i08)) (coin::spearman_test and coin::wilcoxon_test, respectively). By using the permutation version, we do not have to make assumptions on the underlying distributions, hence making the tests more robust.

**Important: This is still a work in progress. Failures can be expected. More functions will be added in the future. Furthermore, I have to add more unit tests to test robustness of some functions**. 


To do: 

- Extend this readme file and description.
- Create more informative error message when data has <2 unique values in association study
- Create function that can test multiple pairs of associations (maybe merge with association_study?)
- Merge the function association_study with association_study_long?
- Improve random forest function, implement unit tests
- Create package vignettes

