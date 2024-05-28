### Get start and end times of windows of adjacent significant samples
### Fall 2021
### adam.milton.morgan@gmail.com

# Given a time series (class: vector) of 1s (significant) and 0s (not significant),
# return either a list or data.frame ("output.class") summarizing the significant windows --
# i.e., temporally adjacent samples that are all significant.

get.significant.windows <- 
  function(sig.vals,
           .times = NULL,
           .samples = NULL,
           .sample.labels = NULL,
           .alternative.labels = NULL,
           output.class = c('list','data.frame')[1],
           include.duration = FALSE,
           .exclude.times.before.ms = NULL,
           .exclude.times.after.ms = NULL,
           .exclude.sig.durations.under.ms = NULL,
           .output.units = c('times','samples','alternative.labels')[1]){
    
    library('dplyr')
    # Function to convert between samples, sample labels, and times
    source(paste0(path,'/analysis/R/functions/time_convert.R'))
    
    # Create ordinal indices if not given
    if(sum(as.numeric(is.null(.times)), 
           as.numeric(is.null(.samples)), 
           as.numeric(is.null(.alternative.labels)), 
           as.numeric(is.null(.sample.labels))) != 3){
      message('WARNING: One (and only one) of these must be given: .times, .samples, .sample.labels, .alternative.labels. Assigning samples starting at 1.')
      .samples = 1:length(sig.vals)
    }
    if(!is.null(.times)){
      .samples <- time.convert(.times, "times", "samples")
    }
    if(!is.null(.sample.labels)){
      .samples <- time.convert(.sample.labels, "sample.labels", "samples")
    }
    if(!is.null(.alternative.labels)){
      .samples <- 1:length(.alternative.labels)
    }
    
    # If length(sig.vals) == 0, return message but set sig.vals to a vector of 0s
    if(length(sig.vals) == 0){
      message("WARNING: NO INPUT DATA!")
    }
    
    # If all values NA
    if(all(is.na(sig.vals))){
      message("WARNING (get.significant.windows): All input vals are NA. Treating as 0s")
      sig.vals <- rep(0, times = length(sig.vals))
    }
    
    # Exclude times before .exclude.times.before.ms
    if(!is.null(.exclude.times.before.ms)){
      exclude.these <- which(.samples < time.convert(.exclude.times.before.ms, "times", "samples"))
      if(length(exclude.these) > 0){
        sig.vals <- sig.vals[-exclude.these]
        .samples <- .samples[-exclude.these]
      }
      rm(exclude.these)
    }
    
    # Exclude times after .exclude.times.after.ms
    if(!is.null(.exclude.times.after.ms)){
      exclude.these <- which(.samples > time.convert(.exclude.times.after.ms, "times", "samples"))
      if(length(exclude.these) > 0){
        sig.vals <- sig.vals[-exclude.these]
        .samples <- .samples[-exclude.these]
      }
      rm(exclude.these)
    }
    
    # Error message
    if(length(sig.vals) == 0){message("Error: No sig values to analyze (probably because .exclude.times.before.ms is positive instead of negative).")}
    
    # Create list of groups of adjacent significant samples
    sig.windows <- list()
    window.number <- 0
    # Manually code first entry if it starts on first row
    if(sig.vals[1] == 1){
      window.number <- 1
      sig.windows[['w1']] <- c(.samples[1], NA)
    }
    # Loop thru rows and add list entries for each new window of adjacent sig samples
    for(sample.loop in 2:length(sig.vals)){ # sample.loop = 769
      if((sig.vals[sample.loop] == 1) & 
         (sample.loop != length(sig.vals))){ # if just the last value is a lone sig window it'll mess things up, just ignore
        if(sig.vals[sample.loop - 1] == 0){
          # Initialize entry
          window.number <- window.number + 1
          sig.windows[[paste0('w',window.number)]] <- 
            c(.samples[sample.loop], NA)
        }
      }else{ # i.e., if current row is not significant
        if(sig.vals[sample.loop - 1] == 1){
          sig.windows[[paste0('w',window.number)]][2] <- 
            .samples[sample.loop]
          names(sig.windows[[paste0('w',window.number)]]) <- c('start.time','end.time')
        }
      }
      # If the last row is still in a sig window, add last time
      if((sample.loop == length(sig.vals)) & 
         sig.vals[sample.loop] == 1 &
         sig.vals[sample.loop - 1] == 1){
        sig.windows[[paste0('w',window.number)]][2] <- 
          .samples[sample.loop]
        names(sig.windows[[paste0('w',window.number)]]) <- c('start.time','end.time')
      }
    }
    
    # Calculate duration if .exclude.sig.durations.under.ms == TRUE (N.B., may override "include.duration = FALSE")
    if(!is.null(.exclude.sig.durations.under.ms)){
      include.duration <- TRUE
    }
    
    if(include.duration & (length(sig.windows) > 0)){
      for(i in length(sig.windows):1){
        sig.windows[[i]]['duration'] <- 
          sig.windows[[i]]['end.time'] - sig.windows[[i]]['start.time']
        # Remove if below .exclude.sig.durations.under.ms
        if(!is.null(.exclude.sig.durations.under.ms)){
          if(sig.windows[[i]]['duration'] < time.convert(.exclude.sig.durations.under.ms, "times", "samples")){
            sig.windows[[i]] <- NULL
          } # if(sig.windows[[i]]['duration'] < .exclude.sig.durations.under.ms)
        } # if(!is.null(.exclude.sig.durations.under.ms))
      }; rm(i)
    } # if(include.duration & (length(sig.windows) > 0))
    
    # Convert output to specified units
    if(length(sig.windows) > 0){
      # If given in .alternative.labels (e.g., frequencies for power spectral density data), return to those labels
      if(!is.null(.alternative.labels)){
        for(i in 1:length(sig.windows)){
          current.names <- names(sig.windows[[i]])
          for(j in 1:length(sig.windows[[i]])){
            sig.windows[[i]][j] <- .alternative.labels[sig.windows[[i]][j]]
          }; rm(j)
          names(sig.windows[[i]]) <- current.names
        }; rm(i, current.names)  
      }else{ # if(!is.null(.alternative.labels)){
        if(.output.units != 'samples'){
          for(i in 1:length(sig.windows)){
            current.names <- names(sig.windows[[i]])
            sig.windows[[i]] <- time.convert(sig.windows[[i]], "samples", .output.units)
            names(sig.windows[[i]]) <- current.names
          }; rm(i, current.names)
        } # if(.output.units != 'samples'){
      } # if(!is.null(.alternative.labels)){}else{
    } # if(length(sig.windows) > 0){
    
    # Convert output to specified class
    if(output.class == 'data.frame'){
      sig.windows <- data.frame(bind_rows(sig.windows))
      if(nrow(sig.windows) == 0){
        sig.windows <- data.frame('start.time' = NA, 'end.time' = NA, 'duration' = NA)[0,]
      }
    }
    
    # End
    return(sig.windows)
  }