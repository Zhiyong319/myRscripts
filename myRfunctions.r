#
# Collect all the useful functions
#
# Author: Zhiyong Wu (zhiyong319@gmail.com)

match_rows <- function(row_match,row_raw) {
  indicies <- vector(length=length(row_raw))
  for (i in 1:length(row_raw)) {
    indicies[i] <- which(row_match==row_raw[i])
  }
                     
  return(indicies)
}