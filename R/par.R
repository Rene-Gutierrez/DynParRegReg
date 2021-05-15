#' Smallest Partition
#'
#' Creates a new partition with the smallest number of elements based on a
#' partition and a maximum number of elements per part.
#'
#' @param P A list with a partition.
#' @param M Maximum Number of Components.
#'
#' @return A list with the new partition.
#'
#' @author Rene Gutierrez Marquez

#' @export


par <- function(P, M){
  # Gets the length of the Connected Components
  n <- length(P)
  # Creates a New Partition
  nP <- list()
  # Goes through every connected Component Randomly
  # ranOrd <- sample(x = n)
  j      <- 1
  for(i in 1:n){
    # Checks if there is a new Component
    if(j > length(nP)){
      # Fills the Component
      nP[[j]] <- P[[i]]
    # Checks if the Component Can fit more elements
    } else if(length(nP[[j]]) + length(P[[i]]) <= M) {
      # Adds the elements
      nP[[j]] <- c(nP[[j]], P[[i]])
    } else {
      # Creates a New Component
      j <- j + 1
      nP[[j]] <- P[[i]]
    }
  }
  # Returns the New Partition
  return(nP)
}
