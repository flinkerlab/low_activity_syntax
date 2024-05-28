### Return the minimum and the maximum of a set of values (like "range()" but with labels and bettern na.rm control)
### Summer 2023
### adam.milton.morgan@gmail.com

min.max <- function(vals,
        .na.rm = FALSE){
  # Get original class
  og.class <- class(vals)
  # Make sure it's a list
  if(og.class != "list"){
    vals <- list(vals)
  }
  
  # Loop thru list entries and return vector c('max','min')
  out <- list()
  if(is.null(names(vals))){
    loop.thru.these <- 1:length(vals)
  }else{
    loop.thru.these <- names(vals)
  }
  
  # Names of values
  for(i in loop.thru.these){
    out[[i]] <- c('min' = min(vals[[i]], na.rm = .na.rm), 
                  'max' = max(vals[[i]], na.rm = .na.rm))
  }
  
  # Return to og class
  if(length(out) == 1 & og.class != "list"){
    out <- out[[1]]
  }
  
  # Return
  return(out)
} # max.min()

