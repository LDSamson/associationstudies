# During my PhD I wrote several functions with the task to loop over variables in a dataset.
# these functions were only slightly different from each other; therefore
# they are difficult to maintain.
# some of these functions are: perform.single.pair.test ,
# all.pairs.to.test, associations.per.ID
# This script below was used to test several different ways to loop over data, and
# to check whether the number of looping functions can be reduced effectively
# so that the package is better to maintain. The function test_association
# shold become the core function of the package.

# I thus check whether pairs of associations can also be tested effectively with
# the test_association function

# load data and functions from chapter 6 for these tests
data.to.analyze <- biomarker.slope.auc.data %>%
  select(ID, Batch, Biomarker, AUC.scaled) %>%
  pivot_wider(names_from = Biomarker, values_from = AUC.scaled)
all.pairs.to.test <- variable_pair_combinations(data.to.analyze, -c(ID, Batch))

# old version of loop
b <- apply(
  all.pairs.to.test,
  MARGIN = 1,
  perform.single.pair.test,
  data.set = data.to.analyze,
  stratum = c("Batch"),
  resampling.number = 10^3
) %>%
  do.call("rbind", .) %>%
  perform.BH.selection()

# new version, base r
base_r_version <- apply(all.pairs.to.test,
                        MARGIN = 1,
                        FUN = function(x){
                          test_association(
                            dataset = data.to.analyze,
                            response.var = x[1],
                            explanatory.var = x[2],
                            stratum = c("Batch"),
                            n.resample = 10^3
                          )
                        }) %>%
  do.call("rbind", .)

# new version, pmap_dfr:
all.pairs.to.test <- variable_pair_combinations(data.to.analyze, -c(ID, Batch))
pairs <- as.data.frame(all.pairs.to.test)
a <- function(){pmap_dfr(pairs,
                         function(V1, V2) {test_association(
                           dataset = data.to.analyze,
                           response.var = V1,
                           explanatory.var = V2,
                           stratum = c("Batch"),
                           n.resample = 10^3
                         )}) }
a <- a %>% arrange(desc(rho))

# benchmark results:
library(microbenchmark)
mbm <- microbenchmark::microbenchmark(
  base = apply(all.pairs.to.test,
               MARGIN = 1,
               FUN = function(x){
                 test_association(
                   dataset = data.to.analyze,
                   response.var = x[1],
                   explanatory.var = x[2],
                   stratum = c("Batch"),
                   n.resample = 10^3
                 )
               }) %>%
    do.call("rbind", .),
  alt = pmap_dfr(pairs,
                 function(V1, V2) {test_association(
                   dataset = data.to.analyze,
                   response.var = V1,
                   explanatory.var = V2,
                   stratum = c("Batch"),
                   n.resample = 10^3
                 )}),
  old = apply(
    all.pairs.to.test,
    MARGIN = 1,
    perform.single.pair.test,
    data.set = data.to.analyze,
    stratum = c("Batch"),
    resampling.number = 10^3
  ) %>%
    do.call("rbind", .),
  times = 20
)
summary(mbm)
#expr      min       lq     mean   median       uq       max neval cld
#1 base 2.868282 3.368535 3.580856 3.515144 3.775343  4.593722    20  a
#2  alt 2.694439 3.391336 3.849051 3.627324 3.868630  7.016552    20  a
#3  old 5.640612 7.071019 8.911306 7.408381 8.947996 21.495821    20   b


# C/ old version is slowest. Both base version and pmap version have similar performance.
# thus, don't use function perform.single.pair.test anymore: use test_association.
#
# all.pairs.to.test should be analyzed and performed


