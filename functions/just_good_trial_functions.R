### Functions for subsetting data to just good trials
### Spring 2023
### adam.milton.morgan@gmail.com

###
### Picture naming
###

just.good.pic.naming <- function(.data,
                                 .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                 .remove.outliers = TRUE,
                                 .remove.over.3000ms = TRUE,
                                 .rt.quantile.range = c(.025,.95),
                                 .rt.range.ms = NULL){
  
  # .data: a dataframe with at least columns for word, substitution_error, pos, production_latency, and case
  # .words: the words to subset to
  # .remove.outliers: exclude trials with RTs greater than or less than specified values (either .rt.quantile.range or .rt.range.ms)
  # .rt.quantile.range: the RT quantiles above and below which to exclude trials; provide either this or the range in RTs (.rt.range.ms; milliseconds)
  # .rt.range.ms: the RTs above and below which to exclude trials; provide either this or the range in quantiles (.rt.quantile.range; proportions)

  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == 'none') & 
                                             is.na(.data$substitution_error) &
                                             (.data$pos == 'n') &
                                             (.data$word %in% .words) &
                                             (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                     (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
}


###
### Sentence data -- first word only
###

just.first.word.sentence.data <- function(.data,
                                          .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                          .remove.over.3000ms = TRUE,
                                          .remove.outliers = FALSE,
                                          .rt.quantile.range = c(.025,.95),
                                          .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == 's') & 
                                   is.na(.data$substitution_error) &
                                   (((.data$dp_or_np == 'np') & (.data$pos == 'n') & (.data$word %in% .words)) | 
                                      ((.data$dp_or_np == 'dp') & (.data$pos == 'd') & (.data$word == 'the'))) &
                                   (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.first.word.sentence.data()


###
### List data -- first word only
###

just.first.word.list.data <- function(.data,
                                      .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                      .remove.over.3000ms = TRUE,
                                      .remove.outliers = FALSE,
                                      .rt.quantile.range = c(.025,.95),
                                      .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == '1') & 
                                   is.na(.data$substitution_error) &
                                   (.data$pos == 'n') &
                                   (.data$word %in% .words) &
                                   (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.first.word.list.data()

