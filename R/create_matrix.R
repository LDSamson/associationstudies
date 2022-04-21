#' Create corrplot matrix
#' function to create FDR or rho matrix for input in the corrplot function of
#' the package corrplot
#' Expects input from p values/FDR values or rhos, in format as table output from
#' association studies (see function association_study) (correlation results)
#'
#' @param data input data (data frame or matrix)
#' @param x first variable
#' @param y second variable
#' @param values values to use
#' @param markers.to.test all biomarkers to test
#'
#' @return matrix for input in corrplot function
#' @export
#'
create_matrix <- function(data, x = "first.var",
                          y = "second.var",values = "rho", markers.to.test){
  matrix.1 <- dplyr::select(data, x, y, values)
  matrix.2 <- dplyr::select(data, y, x, values)
  names(matrix.2) <- names(matrix.2)[c(2,1,3)]
  matrix.full <- dplyr::bind_rows(matrix.1, matrix.2) %>%
    tidyr::pivot_wider(names_from = y, values_from = values)
  matrix.full[[x]] <- factor(matrix.full[[x]], levels = markers.to.test)
  matrix.full <- matrix.full %>%
    dplyr::arrange(!!rlang::sym(x)) %>%
    as.data.frame()
  if(values == "rho"){
    matrix.full <- matrix.full %>%
      # for the middle row (cor=1), replace NA (not tested) with cor=1
      dplyr::mutate_all(.funs = ~ifelse(is.na(.), 1, .))
  }
  if(values == "FDR_selection"){
    matrix.full <- matrix.full %>%
      # for the middle row, replace NA (not tested) with FDR_selection=0
      dplyr::mutate_all(.funs = ~ifelse(is.na(.), 0, .)) %>%
      # ! below will change the FDR selection so that it can be used directly
      # in the corrplot function as if they are p values
      dplyr::mutate_all(.funs = ~ifelse(. == 1, 0.01, 0.99))
  }
  rownames(matrix.full) <- markers.to.test
  matrix.full <- matrix.full %>%
    dplyr::select(markers.to.test) %>%
    as.matrix()
  return(matrix.full)
}
