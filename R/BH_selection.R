## 08-01-2020 LS:
##    - Fixed bug: Now returns warning and a valid table
##      when no associations are found

#' BH selection
#'
#' Function to create a selection based on a set False Discovery Rate
#' I used this function in combination with the association_study function.
#'
#' @param data a data frame
#' @param p_value character value of column containing the p values
#' @param FDR_cutoff numerical value. FDR cutoff
#' @param hist logical. Whether or not to create a histogram
#'
#' @return data frame with added columns with FDR and FDR selection
#' @importFrom rlang .data
#' @export
BH_selection <- function(data, p_value = "p.value",
                         FDR_cutoff = 0.15, hist=FALSE){
  # calculate FDR
  n.of.tests <- nrow(data)
  data.sorted <- data %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::arrange(!!rlang::sym(p_value)) %>%
    dplyr::mutate(
      'order.of.tests' = dplyr::row_number(),
      'FDR'            = n.of.tests*!!rlang::sym(p_value)/.data$order.of.tests,
      FDR              = round(.data$FDR, digits = 4)
    )
  if(min(data.sorted$FDR) > FDR_cutoff){
    warning(paste0("No associations detected with FDR of ", FDR_cutoff))
    data.with.BH.selection <- data.sorted %>%
      dplyr::mutate(`FDR_selection` = 0) %>%
      dplyr::select(-.data$order.of.tests)
  }else{
    # Define max number of detected associations, based on FDR cutoff
    max.n.of.tests <- max(data.sorted[data.sorted$FDR <= FDR_cutoff,
                                      "order.of.tests"])
    data.with.BH.selection <- data.sorted %>%
      dplyr::mutate(
        FDR_selection = ifelse(.data$order.of.tests <= max.n.of.tests, 1, 0)
      ) %>% dplyr::select(-.data$order.of.tests)
  }
  # (optional): print p value histogram
  if(isTRUE(hist)){plot_p_histogram(data.with.BH.selection)}
  return(data.with.BH.selection)
}
