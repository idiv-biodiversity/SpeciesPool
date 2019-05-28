#' Rarefaction curves
#'
#' Builds an empirical rarefaction curve
#' @author Francesco Maria Sabatini
#' @author Helge Bruelheide
#' @param x A vector of species richness
#' @param areas A vector of plot sizes
#' @param n numeric - number of intervals of rarefaction curve
#' @return Returns a dataframe of cumulative richness and pooled plot areas
#' @export



###### SAR for type IIIb
SAR.IIIb <- function(x, areas, n=10){
  ll <- length(areas)
  ks <- c(unique(round(10^(seq(log10(1),log10(ll-1), length=10)))),ll)
  choo <- choose(ll,ks)
  choo2 <- choo
  choo2[choo>n] <- n
  output <- list()
  #mean.rich <- rep(NA,length(ks))
  rich.out <- matrix(NA, nrow=length(ks), ncol=n, dimnames=list(ks,1:n))
  areas.out <- matrix(NA, nrow=length(ks), ncol=n, dimnames=list(ks,1:n))
  for(i in 1:length(ks)){
    if( (choo >= n)[i]){
      tmp <- replicate(n, sample(ll,ks[i],replace = F))
      if(is.vector(tmp)) tmp <- matrix(tmp, nrow=1)
    } else {
      tmp <- combn(ll, ks[i])
    }

    #rich.j <- rep(NA,choo2[i] )
    for(j in 1:choo2[i]){
      rich.out[i,j] <- sum(colSums(x[tmp[,j],,drop=F])>0)
      areas.out[i,j] <- sum(areas[tmp[,j]])
    }
    #mean.rich[i] <- mean(rich.j)
  }
  return(data.frame(melt(rich.out, varnames = c("num.pooled.plots", "replicate"), value.name="Pooled.Richness", na.rm = T),
                    Sum.area=melt(areas.out, varnames = c("num.pooled.plots", "replicate"), value.name="Sum.area", na.rm = T)$Sum.area))
}
