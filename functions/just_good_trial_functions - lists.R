### Functions for subsetting data to just good trials
### Spring 2023
### adam.milton.morgan@gmail.com

###
### Picture naming
###

just.good.pic.naming <- function(.data,
                                 .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                 .remove.outliers = TRUE,
                                 .rt.quantile.range = c(.025,.95)){
  # Get input class
  input.class <- class(.data)
  
  # If just a data.frame (e.g., just trial info), make it a list
  if(input.class == 'data.frame'){.data <- list(.data)}
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- lapply(.data,
                  function(x){droplevels(x[(x$case == 'none') & 
                                             is.na(x$substitution_error) &
                                             (x$pos == 'n') &
                                             (x$production_latency < time.convert(3000, 'times', 'samples')) &
                                             (x$word %in% .words) &
                                             (! is.na(x$case)),])})
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data[[1]]$production_latency
    
    # Get cutoff RTs (in samples)
    rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
    rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    
    # Subset to just trials within cutoff RTs
    .data <- lapply(.data,
                    function(x){
                      droplevels(x[(x$production_latency >= rts.samples.lower.bound) &
                                     (x$production_latency <= rts.samples.upper.bound),])
                    })
  } # if remove outliers
  
  # Restore input class
  if(input.class == 'data.frame'){.data <- .data[[1]]}
  
  # Return trimmed data
  return(.data)
}


###
### Sentence data -- first word only
###

just.first.word.sentence.data <- function(.data,
                                          .words = c('chicken','dog','dracula','frankenstein','ninja','nurse')){
  # Get input class
  input.class <- class(.data)
  
  # If just a data.frame (e.g., just trial info), make it a list
  if(input.class == 'data.frame'){.data <- list(.data)}
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- lapply(.data,
                  function(x){
                    droplevels(x[(x$case == 's') & 
                                   is.na(x$substitution_error) &
                                   (((x$dp_or_np == 'np') & (x$pos == 'n') & (x$word %in% .words)) | 
                                      ((x$dp_or_np == 'dp') & (x$pos == 'd') & (x$word == 'the'))) &
                                   (! is.na(x$case)),])})
  
  # Restore input class
  if(input.class == 'data.frame'){.data <- .data[[1]]}
  
  # Return trimmed data
  return(.data)
} # just.first.word.sentence.data()


###
### List data -- first word only
###

just.first.word.list.data <- function(.data,
                                      .words = c('chicken','dog','dracula','frankenstein','ninja','nurse')){
  # Get input class
  input.class <- class(.data)
  
  # If just a data.frame (e.g., just trial info), make it a list
  if(input.class == 'data.frame'){.data <- list(.data)}
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- lapply(.data,
                  function(x){
                    droplevels(x[(x$case == '1') & 
                                   is.na(x$substitution_error) &
                                   (x$pos == 'n') &
                                   (x$word %in% .words) &
                                   (! is.na(x$case)),])})
  
  # Restore input class
  if(input.class == 'data.frame'){.data <- .data[[1]]}
  
  # Return trimmed data
  return(.data)
} # just.first.word.list.data()

