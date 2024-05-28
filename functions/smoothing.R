### Boxcar smoothing function
### Winter 2021
### adam.milton.morgan@gmail.com

smoothing <- function(time.series,
                      sample.labels=NA,
                      n.samples.pre=30, 
                      n.samples.post=n.samples.pre, 
                      na.pad=TRUE){
  
  # Set up labels
  if(all(is.na(sample.labels))){
    labels = names(time.series)
  }else{
    labels = sample.labels
  }
  
  # Convert
  time.series = as.numeric(time.series)
  output = c()
  for(i in (n.samples.pre+1):(length(time.series)-n.samples.post)){
    output = c(output, mean(time.series[(i - n.samples.pre):(i + n.samples.post)],
                            na.rm=TRUE))
  }
  if(na.pad){
    output = c(rep(NA, n.samples.pre), output, rep(NA, n.samples.post))
    # Add labels
    names(output) = labels
  }else{
    names(output) = labels[(n.samples.pre+1):(length(time.series)-n.samples.post)]
  }
  
  # End
  return(output)
}