#' Pairwise likelihood of species co-occurrence
#'
#' Intermediate step to calculate Beals' smoothing
#' @author Francesco Maria Sabatini
#' @author Helge Bruelheide
#' @param dt A data.frame with three columns: releve ID, species ID, and abundance
#' @return Returns a square matrix of the likelihood of pairwise species co-occurrence
#' @export

Mij.calc <- function(dt, verbose=F) {
  if(ncol(dt) !=3) stop("dt should have three columns:  releve, species,abundance")
  colnames(dt) <- c( "RELEVE", "SPECIES","presabs")
  if(!is.factor(dt$SPECIES)) {
    if(verbose) print("converting species column to factor")
    dt$SPECIES <- factor(dt$SPECIES)
  }
  if(!is.factor(dt$RELEVE)) {
    if(verbose) print("converting releve column to factor")
    dt$RELEVE <- factor(dt$RELEVE)}
  if(!is.numeric(dt$presabs)) {stop("The column abundance should be numeric")}
  if(!any(dt$presabs %in% c(0,1))) {
    if(verbose) print("converting abundances to presence\absence")
    dt$presabs <- as.numeric(dt$presabs>0*1)}
  if(verbose) print("Calculating matrix Mij")
  DT.matrix.pa <- Matrix::sparseMatrix(as.integer(dt$RELEVE),
                               as.integer(dt$SPECIES),
                               x = dt$presabs,
                               dimnames=list(levels(dt), levels(dt$SPECIES)))
  Mij <- Matrix::crossprod(DT.matrix.pa, DT.matrix.pa)
  #Mij <- sweep(Mij, 2, diag(Mij), FUN="/")
  Mij <- sweep(Mij, 2, Mij[seq(1, length(Mij), by=ncol(Mij)+1)], FUN="/") #alternative to above, to subset long integers
  return(Mij)
}
