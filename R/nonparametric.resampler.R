nonparametric.resampler <- function(data){
  resampled.data <- dplyr::sample_n(data, size = nrow(data))
  return(resampled.data)
}