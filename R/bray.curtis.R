#' Bray-Curtis dissimilarity
#'
#' Calculates Bray-Curtis dissimilarity between a target plot, and all the other relevés in the dataset
#' @author Francesco Maria Sabatini
#' @author Helge Bruelheide
#' @param input.data A species x relevés matrix with abundance values
#' @param target.plot.ID ID or label of target plot
#' @return A vector of Bray-Curtis dissimilarities between the target and all the other relevés in input.data
#' @export


bray.curtis <- function(input.data, target.plot.ID){
  # input.data is the matrix with beals values of all plots
  # including the target plot
  # target.plot.ID is the RELEVE_NR of the target plot
  # Bray-Curtis dissimilarity:
  # d[jk] = (sum abs(x[ij]-x[ik]))/(sum (x[ij]+x[ik]))
  target.plot <- input.data[dimnames(input.data)[[1]]==target.plot.ID,]
  #input.data <- input.data[dimnames(input.data)[[1]]!=target.plot.ID,]
  input.data.minus <- sweep(input.data, 2, target.plot, FUN = "-")
  input.data.minus <- abs(input.data.minus)
  input.data.plus <- sweep(input.data, 2, target.plot, FUN = "+")
  bray <- rowSums(input.data.minus)/rowSums(input.data.plus)
  return(bray)
}
