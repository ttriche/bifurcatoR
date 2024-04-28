#' Our own initz with fuzzy c means if kmeans or hclust fails in mixR::initz
#' 
#' This was Zach Madaj's idea and it fixes the underlying issue in mixR, which 
#' as of April 2024, remained in spite of previous attempts to patch it. Turns
#' out this is the least worst way to deal with a recurrent issue in 1D models.
#'
#' nasty kludge:
#' # unlockBinding("initz", as.environment("package:mixR"))
#' # assign("initz", initz2, "package:mixR")
#'
#' Hopefully we can eliminate the nasty kludge by contributing this to mixR.
#'
#' @param x             the data to fit
#' @param ncomp         how many components to model 
#' @param init.method   initialization method ('fuzzy' is new to initz2)
#' @param minprop       minimum proportion of data in a component before fuzzing
#'
#' @return              a list with values "pi", "mu", and "sd" of length ncomp
#'
#' @importFrom ppclust  fpppcm
#' 
#' @seealso ppclust::fpppcm
#' @seealso mixR::initz
#'
#' @export
#'
initz2 <-  function(x, ncomp, init.method = c("kmeans", "hclust","fuzzy"), minprop=0.05) {
  
  init.method <- match.arg(init.method)

  # check if 'x' is a matrix (from grouped data)
  if(is.matrix(x)) {
    x <- reinstate(x)
  }

  if(init.method == "kmeans") {
    a <- kmeans(x, centers = ncomp,nstart = 1)$cluster
    if(any(table(a) < (length(x) * minprop))){
      a <- fpppcm(x,centers = ncomp)$cluster
    }
  } else {
    a <- cutree(hclust(dist(x)), ncomp)
    # a <- fpppcm(x,centers = ncomp)$cluster
  } 

  res <- list()
  for(i in 1:ncomp) {
    res[[i]] <- x[a == i]
  }
  count <- sapply(res, length)
  pi <- count / sum(count)
  mu <- sapply(res, mean)
  sd <- sapply(res, sd)
  
  order <- order(mu)
  
  pi <- pi[order]
  mu <- mu[order]
  sd <- sd[order]
  ret <- list(pi = pi, mu = mu, sd = sd)
  return(ret)

}
