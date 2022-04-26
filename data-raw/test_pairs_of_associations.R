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



#### tryout scientific notation function below:

#' Scientific notation
#'
#' @param x input (vector of numerical values)
#' @param n.digits number of digits to round the results
#' @param scientific.above numerical value, above which the scientific notation hsold be used
#'
#' @return
#'
#' @examples
scient_func <- function(x, n.digits = 3, scientific.above = 1000){
  if(is.na(n.digits)) n.digits <- 0
  if(n.digits<0){
    stop("Positive integer value needed for n.digits")}
  if(scientific.above<0){
    stop("Positive integer value needed for scientific.above")}

  xabs <- format(abs(x), scientific = FALSE)
  return.original <- abs(x) < scientific.above & abs(x) > 1/scientific.above
  xchar <- trimws(as.character(xabs))
  xchar <- gsub("\\.", "", xchar)
  xchar <- gsub("^0+", "", xchar)

  base <- substr(xchar, 1, 1)
  digits <- substr(xchar, 2, 2 + n.digits) # one more digit for rounding later
  base.num <- formatC(as.numeric(paste(base, digits, sep = ".")), digits = n.digits, format = "f")
  return.vals <- paste0(base.num, "*10^", floor(log10(abs(x))))
  return.vals[return.original] <- formatC(abs(x[return.original]), digits = n.digits, format = "f")
  return.vals[x<0] <- paste0("-", return.vals[x<0])
  return.vals
}
#x <- c(20000, 15000, 25000)
x <- c(20.15389, 2.1, 300, 21678, 2010, -1000, -1.23*10^4, 0.000001234, 0.167)
scient_func(x, n.digits = 2, scientific.above = 1000)
sumval_table_format(abs(x))

# other option would be to directly use outcomes of formatC function, change the
a <- formatC(x, digits = 2, format = "g")
a <- gsub("e\\+", "*10^", a)
a <- gsub("e\\-", "*10^-", a)
a <- gsub("\\^0", "\\^", a)
n.digits <- 2

a <- scient_func(c(20000, 15000, 25000), n.digits = 2, scientific.above = 1000)
exponent <- stringr::str_extract(a[1], "\\*.*$")
if(is.na(exponent)) exponent <- ""
a <- gsub("\\*.*$", "", a)
paste0(a[1], "(", a[2], "-", a[3], ")", exponent)
a <- c(rnorm(100, 10000, 500), 0, 0)
sumval_table_format(a, n.digits = 2, conf.level = 0.95)
sumval_table_format(rnorm(100, 10, 4))
a <- abs(rnorm(100, 1, 0.4))
gm_mean(a)
sumval_table_format(a*1000, n.digits = 2)
formatC(x, digits = 3, format = "g", drop0trailing = TRUE)

## scientific simple
scient_func_simple <- function(x, n.digits = 3){
  if(is.na(n.digits)) n.digits <- 0
  if(n.digits<0){
    stop("Positive integer value needed for n.digits")}
  xchar <- formatC(x, format = "g", digits = n.digits, drop0trailing = T)
  xchar <- gsub("e\\+|e\\+0", "*10^", xchar)
  xchar <- gsub("e\\-|e\\-0", "*10^-", xchar)
  xchar
}

# you can combine exponent here, but edge case problem is when exponent
# in confidence interval is different from the geomean value exponent!
x <- rnorm(100, 1000000, 1000000)
x <- x[x > 0]
x <- c(20.15389, 2.1, 300, 21678, 2010, -1000, -1.23*10^4, 0.000001234, 0.167)
x <- abs(x)

gm_table_format <- function(x, n.digits = 2){
  xchar <- scient_func_simple(gm_mean(x, na.rm = T, conf.level = 0.95),
                              n.digits = n.digits)
  exponent <- stringr::str_extract(xchar, "\\*.*$")
  exponent[is.na(exponent)] <- ""
  xchar <- gsub("\\*.*$", "", xchar)
  if(length(unique(exponent)) != 1){
    exp_val <- as.numeric(gsub("\\*10\\^", "", exponent))
    exp_val[is.na(exp_val)] <- 0
    exp_value_diff <- exp_val - exp_val[1]
    xchar <- formatC(as.numeric(xchar)*10^exp_value_diff, format = "f", digits = n.digits)
    }
  return.vals <- paste0(xchar[1], "(", xchar[2], "-", xchar[3], ")", exponent[1])
  return.vals
}
gm_table_format(abs(x), n.digits = 4)
sumval_table_format(abs(x), n.digits = 3)

# problematic edge case:
#a <- abs(rnorm(100, 10000, 10000))
#b <- scient_func_simple(gm_mean(a, conf.level = 0.95))
