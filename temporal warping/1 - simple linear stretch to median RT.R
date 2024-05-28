### Linearly interpolate the middle section (150ms post stim to 150ms pre speech) of all trials to have the same RT (i.e., median task RT)
### May 2023
### adam.milton.morgan@gmail.com

###
### Readme
###


###
### Setup
###

### Packages
library('dplyr') # for data organization, including bind_rows()

# Clean up
rm(list=ls())
cat("\014")
message("Begin stretching trials. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 4
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 5
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 15
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 40
}

### Lemmas
# Function to convert between samples, sample labels, and times
source(paste0(path,'/analysis/R/functions/time_convert.R'))
# Function for rolling average smoothing
source(paste0(path,'/analysis/R/functions/smoothing.R'))
# Close any old parallel backends
source(paste0(path,'/analysis/R/functions/unregister_dopar.R'))
# Subset to just good picture naming trials
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add colored text to plots
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))

set.seed(seed = 404)

# Data type:
for(freq.band.loop in c('high_gamma','beta')){
  # freq.band.loop = c('high_gamma','beta')[1]
  freq.band.read.path <- list('high_gamma' = 'z_scored', 'beta' = 'z_scored_beta')[[freq.band.loop]]
  
  # Type of hilbert transform
  hilbert.type <- 'multiband'
  
  ## Metadata
  patients = c('Patient001',
               'Patient002',
               'Patient003', # slow PN RTs; check closely
               'Patient004',
               'Patient005',
               'Patient006', # slow PN RTs; check closely
               'Patient007',
               'Patient008', # slow PN RTs and generally weird
               'Patient009',
               'Patient010')
  
  ## Tasks
  tasks <- c('pn','sp','lp')
  
  ## Plot stuff:
  load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))
  text.size.big = 1.6
  text.size.med = 1.4
  text.size.small = 1.2
  zoom = 2
      
  # Elec MNI coords, etc.
  elec.info <- 
    read.csv(paste0(path,
                    '/analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
  elec.info$use.these.elecs <- 
    as.numeric((! elec.info$region_clinical %in% c('Unknown',
                                                   '',
                                                   'NaN',
                                                   'Left-Cerebral-White-Matter',
                                                   'Left-Inf-Lat-Vent',
                                                   'Right-Cerebral-White-Matter')) &
                 (elec.info$bad_elec == 0) &
                 (elec.info$visual_elec == 0) &
                 (elec.info$active == 1))
  rownames(elec.info) <- elec.info$patient_elec
  use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec
  
  
  ## Get all trial info across patients in one dataframe to analyze RTs:
  # Define range of RTs per patient to use
  rt.quantile.range <- c(.025, .95) # .025 just to get rid of weird fast trials (e.g., second of two 'nurse' trials in a row?) and .8 because long tail, probably reflecing lots of different things, so get rid of a bunch. shouldn't matter long-run, still a decent amount of data for NMF clustering
  
  ### Load trial labels
  # Initialize storage
  trial.info <- list()
  for(task.loop in c(tasks,'all')){
    trial.info[[task.loop]] <- list()
  }; rm(task.loop)
  
  # Loop thru patients and load data
  for(patient in rev(patients)){
    # patient = patients[1]
    print(patient)
    
    # Load data
    trial.info[['all']][[patient]] <- read.csv(paste0(path,
                                                      'data/',patient,'/data/epoched/',patient,'_trial_labels_without_bad_trials.csv'))
    trial.info[['all']][[patient]]$patient <- patient
    trial.info[['all']][[patient]]$rt <- time.convert(trial.info[['all']][[patient]]$production_latency, 
                                                      "samples", "times")
    
    # Separate by task
    trial.info[['pn']][[patient]] <- 
      just.good.pic.naming(trial.info[['all']][[patient]], 
                           .remove.outliers = TRUE, 
                           .rt.quantile.range = rt.quantile.range)
    trial.info[['sp']][[patient]] <- 
      just.first.word.sentence.data(trial.info[['all']][[patient]], 
                                    .remove.outliers = TRUE, 
                                    .rt.quantile.range = rt.quantile.range)
    trial.info[['lp']][[patient]] <- 
      just.first.word.list.data(trial.info[['all']][[patient]], 
                                .remove.outliers = TRUE, 
                                .rt.quantile.range = rt.quantile.range)
  }; rm(patient)
  
  # Clean up
  trial.info$all <- NULL
  trial.info[['pn']] <- data.frame(bind_rows(trial.info[['pn']]))
  trial.info[['sp']] <- data.frame(bind_rows(trial.info[['sp']]))
  trial.info[['lp']] <- data.frame(bind_rows(trial.info[['lp']][patients[patients != 'Patient006']]))
  
  ## Get median RT by tasks -- i.e., target RT for all adjusted trials
  half.median.rt.samples <- round(sapply(trial.info, function(x){median(x$production_latency)}) / 2)
  median.rt.samples <- 2 * half.median.rt.samples # double rounded halves so always even
  
  # Windows of data keep (i.e., not stretch) pre-stim and post-speech
  stimulus.samples.range <- c(time.convert(-500, "times", "samples"), time.convert(150, "times", "samples"))
  production.samples.range <- c(time.convert(-150, "times", "samples"), time.convert(1200, "times", "samples"))
  stimulus.samples <- stimulus.samples.range[1]:stimulus.samples.range[2]
  production.samples <- production.samples.range[1]:production.samples.range[2]
  stimulus.sample.labels <- time.convert(
    stimulus.samples,
    "samples", "sample.labels")
  production.sample.labels <- time.convert(
    production.samples,
    "samples", "sample.labels")
  
  
  
  ### Read in data from both time-locks
  stretch.data <- list()
  patient.trial.info <- list()
  elec.means.pre.post <- list()
  elec.mean.maxes.pre.post <- list()
  for(patient in rev(patients)){
    # patient = patients[1]
    
    # Initialize storage
    data <- list()
    keep.data.stimulus <- list()
    keep.data.production <- list()
    
    for(lock.loop in c('locked_to_stimulus_onset','locked_to_production_onset')){
      #lock.loop <- c('locked_to_production_onset','locked_to_stimulus_onset')[2]
      
      # Storage
      data[[lock.loop]] <- list()
      
      # Progress update
      message(patient,": Loading ",lock.loop," data. ",Sys.time())
      
      # Number of linguistic columns
      n.lx.cols = read.table(paste0(path,'data/',patient,'/data/epoched/',patient,'_n_linguistic_columns.txt'))[1,1]
      
      # Where data?
      data.dir <- paste0(path,'data/',
                         patient, 
                         '/data/epoched/elec_data_without_bad_trials/',hilbert.type,'/',
                         lock.loop,'/', # always train on data locked to stimulus onset since we care about early stages!
                         freq.band.read.path,'/')
      
      # Read in data
      message('Attach...')
      attach(paste0(data.dir, patient,' elec data.RData')) # reads in list "all.data"
      data[[lock.loop]][['all']] <- all.data
      detach()
      message('...detach!')
      
      ## Limit trials to inner range of RTs
      data[[lock.loop]][['pn']] <- lapply(data[[lock.loop]][['all']], function(x){
        just.good.pic.naming(x, .remove.outliers = TRUE, .rt.quantile.range = rt.quantile.range)})
      data[[lock.loop]][['sp']] <- lapply(data[[lock.loop]][['all']], function(x){
        just.first.word.sentence.data(x, .remove.outliers = TRUE, .rt.quantile.range = rt.quantile.range)})
      data[[lock.loop]][['lp']] <- lapply(data[[lock.loop]][['all']], function(x){
        just.first.word.list.data(x, .remove.outliers = TRUE, .rt.quantile.range = rt.quantile.range)})
      
      # Clean up
      data[[lock.loop]][['all']] <- NULL
    }; rm(lock.loop)    
    
    ## Set aside the "keep" data -- pre-stim and post-speech onsets (with a little extra from the middle)
    patient.trial.info[[patient]] <- list()
    for(task.loop in tasks){  
      # task.loop = tasks[1]
      
      # Store the trial info (linguistic columns) for this patient
      patient.trial.info[[patient]][[task.loop]] <-
        data[['locked_to_production_onset']][[task.loop]][[1]][,1:n.lx.cols]
      
      # Get the "keep" data
      keep.data.stimulus[[task.loop]] <- 
        lapply(data[['locked_to_stimulus_onset']][[task.loop]], function(x){
          x <- x[,stimulus.sample.labels]
          colnames(x) <- paste0('stim.',colnames(x)) # make unique column names to avoid duplicates
          return(x)
        })
      keep.data.production[[task.loop]] <-
        lapply(data[['locked_to_production_onset']][[task.loop]], function(x){
          x <- x[,production.sample.labels]
          colnames(x) <- paste0('prod.',colnames(x)) # make unique column names to avoid duplicates
          return(x)
        })
      
      # Remove the linguistic trial data from the data to be stretched
      for(lock.loop in c('locked_to_stimulus_onset','locked_to_production_onset')){
        data[[lock.loop]][[task.loop]] <-
          lapply(data[[lock.loop]][[task.loop]], function(x){
            x[,-c(1:n.lx.cols)]
          })  
      }; rm(lock.loop)
      
      # Remove the "keep" data (and any even earlier/later than stimulus/production) from the data to be stretched
      data[['locked_to_stimulus_onset']][[task.loop]] <- 
        lapply(data[['locked_to_stimulus_onset']][[task.loop]], function(x){
          x[,-which(time.convert(colnames(x), "sample.labels", "samples") <= max(stimulus.samples))]
        })
      data[['locked_to_production_onset']][[task.loop]] <- 
        lapply(data[['locked_to_production_onset']][[task.loop]], function(x){
          x[,-which(time.convert(colnames(x), "sample.labels", "samples") >= min(production.samples))]
        })
      
    }; rm(task.loop)
    
    
    ### Set up metadata for stretching data
    # Electrodes
    elecs <- names(data$locked_to_stimulus_onset$pn)
    
    # Target number of samples
    stim.n.samples.out <- half.median.rt.samples - stimulus.samples.range[2]
    prod.n.samples.out <- half.median.rt.samples - (-production.samples.range[1])
    
    # Save names for the final, stretched sample labels 
    stretch.sample.labels <- list()
    for(task.loop in tasks){ # task.loop = tasks[1]
      stretch.sample.labels[[task.loop]] <- list()
      # Get samples
      stretch.sample.labels[[task.loop]][['locked_to_stimulus_onset']] <- c(
        stimulus.samples.range[1]:(median.rt.samples[[task.loop]] + production.samples.range[2] + 1)
      )
      stretch.sample.labels[[task.loop]][['locked_to_production_onset']] <- c(
        stretch.sample.labels[[task.loop]][['locked_to_stimulus_onset']] - median.rt.samples[[task.loop]]
      )
      # Convert to sample labels
      stretch.sample.labels[[task.loop]][['locked_to_stimulus_onset']] <-
        time.convert(stretch.sample.labels[[task.loop]][['locked_to_stimulus_onset']],
                     "samples", "sample.labels")
      stretch.sample.labels[[task.loop]][['locked_to_production_onset']] <-
        time.convert(stretch.sample.labels[[task.loop]][['locked_to_production_onset']],
                     "samples", "sample.labels")
    }; rm(task.loop)
    
    
    ### For each trial, get the "stretch" data as a single list entry -- but *just* the data from 150ms post-stim to 150ms pre-production
    # Take the first half from the stim-locked data and second half from the prod-locked data to cover cases where the RT is longer than the epoch
    stretch.data[[patient]] <- list()
    for(task.loop in tasks){
      # task.loop = tasks[3]
      
      message(patient,': Beginning stretching of ',task.loop,' data. ',Sys.time())
      
      # Loop thru elecs
      stretch.data[[patient]][[task.loop]] <- list()
      for(elec.loop in elecs){
        # elec.loop = elecs[1]
        
        .n.trials <- nrow(patient.trial.info[[patient]][[task.loop]])
        if(.n.trials == 0){
          stretch.data[[patient]][[task.loop]][[elec.loop]]
        }else{
          
          # Loop thru trials and interpolate
          stretch.data[[patient]][[task.loop]][[elec.loop]] <- list()
          for(trial.loop in 1:nrow(patient.trial.info[[patient]][[task.loop]])){
            # trial.loop = 1
            
            # Get all the samples from 1 to production onset
            .rt.samples <- patient.trial.info[[patient]][[task.loop]]$production_latency[trial.loop]
            .all.rt.samples <- (stimulus.samples.range[2] + 1) : (.rt.samples + production.samples.range[1] - 1)
            .stim.locked.samples <- 
              time.convert(.all.rt.samples[c(1:floor(length(.all.rt.samples) / 2))], 
                           "samples", "sample.labels")
            .prod.locked.samples <- 
              time.convert(.all.rt.samples[-c(1:floor(length(.all.rt.samples) / 2))] - .rt.samples,
                           "samples", "sample.labels")
            
            # Verify that the value at t=0 prod-locked is the same as the value at the RT stim-locked
            if(! all(round(data[['locked_to_stimulus_onset']][[task.loop]][[elec.loop]][
              trial.loop, c(.stim.locked.samples[length(.stim.locked.samples)],
                            time.convert(time.convert(.stim.locked.samples[length(.stim.locked.samples)], 
                                                      "sample.labels", "samples") + 1, 
                                         "samples", "sample.labels"))], 6) ==
              round(data[['locked_to_production_onset']][[task.loop]][[elec.loop]][
                trial.loop, c(time.convert(time.convert(.prod.locked.samples[1], 
                                                        "sample.labels", "samples") - 1, 
                                           "samples", "sample.labels"),
                              .prod.locked.samples[1])], 6))){
              message('WARNING!!! Stim- and prod-locked data not the same!')
            }
            
            # Interpolate
            .stim.locked.interpolated.data <-
              approx(x = time.convert(.stim.locked.samples, "sample.labels", "samples"),
                     y = unlist(data[['locked_to_stimulus_onset']][[task.loop]][[elec.loop]][trial.loop, .stim.locked.samples]),
                     n = stim.n.samples.out[[task.loop]])$y
            .prod.locked.interpolated.data <-
              approx(x = time.convert(.prod.locked.samples, "sample.labels", "samples"),
                     y = unlist(data[['locked_to_production_onset']][[task.loop]][[elec.loop]][trial.loop, .prod.locked.samples]),
                     n = prod.n.samples.out[[task.loop]])$y
            
            # Turn into a dataframe
            .stim.locked.interpolated.data <- data.frame(matrix(.stim.locked.interpolated.data, nrow = 1))
            .prod.locked.interpolated.data <- data.frame(matrix(.prod.locked.interpolated.data, nrow = 1))
            colnames(.stim.locked.interpolated.data) <- paste0('stim.', 1:stim.n.samples.out[[task.loop]])
            colnames(.prod.locked.interpolated.data) <- paste0('prod.', 1:prod.n.samples.out[[task.loop]])
            
            stretch.data[[patient]][[task.loop]][[elec.loop]][[trial.loop]] <- 
              cbind(.stim.locked.interpolated.data, .prod.locked.interpolated.data)
            
            rm(.stim.locked.interpolated.data, .prod.locked.interpolated.data, .rt.samples, .all.rt.samples, .stim.locked.samples, .prod.locked.samples)
          }; rm(trial.loop) 
          
          # Recombine trials into dataframe
          stretch.data[[patient]][[task.loop]][[elec.loop]] <- bind_rows(stretch.data[[patient]][[task.loop]][[elec.loop]])
          
          # Add stim-locked prod-locked data back in
          stretch.data[[patient]][[task.loop]][[elec.loop]] <- cbind(
            keep.data.stimulus[[task.loop]][[elec.loop]],
            stretch.data[[patient]][[task.loop]][[elec.loop]],
            keep.data.production[[task.loop]][[elec.loop]]
          )
          colnames(stretch.data[[patient]][[task.loop]][[elec.loop]]) <- stretch.sample.labels[[task.loop]][['locked_to_production_onset']]
        } # if(.n.trials > 0){
        
      }; rm(elec.loop)
    }; rm(task.loop)
    
    
    ### Visual check
    ## Plot random trials
    save.fig.path <- paste0(path,
                            'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',
                            freq.band.loop,
                            '/figures/',patient,'/')
    save.fig.path.trial <- paste0(save.fig.path, 'randomly sampled trials/')
    save.fig.path.mean <- paste0(save.fig.path, 'elec means/')
    dir.create(save.fig.path.trial, showWarnings = FALSE, recursive = TRUE)
    dir.create(save.fig.path.mean, showWarnings = FALSE, recursive = TRUE)
    
    elec.means.pre.post[[patient]] <- list()
    elec.mean.maxes.pre.post[[patient]] <- list()
    for(task.loop in tasks){ 
      # task.loop = tasks[1]
      
      elec.means.pre.post[[patient]][[task.loop]] <- list()
      elec.mean.maxes.pre.post[[patient]][[task.loop]] <- list()
      for(elec.loop in elecs){ 
        # elec.loop = elecs[3]
            
        ## Compare difference between [mean unstretched signal] & [mean stretched signal] as a function of signal: 
        # low signal (noise) the should get lower (noise reduction)
        # high signal should get higher if this stretching approach isn't muddying things
        
        .n.trials <- nrow(patient.trial.info[[patient]][[task.loop]])
        if(.n.trials > 0){ # if they completed any trials this task (Patient006 didn't do listing block)
          
          # Samples to compare
          samples.to.compare <- 
            time.convert((-(median.rt.samples[task.loop] - stimulus.samples.range[2] - 1)):
                           (production.samples.range[1]),
                         "samples", "sample.labels")
          
          # Stretched data
          .stretch.mean <- colMeans(stretch.data[[patient]][[task.loop]][[elec.loop]])
          .stretch.peak.time <- 
            time.convert(names(which.max(abs(.stretch.mean))),
                         "sample.labels", "times")
          
          # Means
          elec.means.pre.post[[patient]][[task.loop]][[elec.loop]] <-
            data.frame('sample.label' = samples.to.compare,
                       'elec' = elec.loop,
                       'task' = task.loop,
                       'localization' = elec.info[elec.loop, 'region_clinical'],
                       'stretch.mean' = abs(.stretch.mean[samples.to.compare]),
                       row.names = NULL)
          
          # Maxes
          elec.mean.maxes.pre.post[[patient]][[task.loop]][[elec.loop]] <- 
            data.frame('patient' = patient,
                       'elec' = elec.loop,
                       'task' = task.loop,
                       'localization' = elec.info[elec.loop, 'region_clinical'],
                       'stretch.max' = max(abs(.stretch.mean[samples.to.compare]), na.rm = TRUE),
                       'stretch.peak.time' = .stretch.peak.time,
                       row.names = NULL)
          
          for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
            # lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')[1]
            
            # Redefine samples to compare if stim-locked
            if(lock.loop == 'locked_to_stimulus_onset'){
              samples.to.compare <- 
                time.convert((stimulus.samples.range[2]):
                               (median.rt.samples[task.loop] + production.samples.range[1] - 1),
                             "samples", "sample.labels")
            } # if(lock.loop == 'locked_to_stimulus_onset'){
            
            # Get mean of OG data locked to lock.loop
            .original.mean <- colMeans(data[[lock.loop]][[task.loop]][[elec.loop]])
            .original.peak.time <- 
              time.convert(names(which.max(abs(.original.mean))),
                           "sample.labels", "times")
            
            # Store original mean and peak
            elec.means.pre.post[[patient]][[task.loop]][[elec.loop]][,paste0('original.mean_',lock.loop)] <-
              abs(.original.mean[samples.to.compare])
            elec.mean.maxes.pre.post[[patient]][[task.loop]][[elec.loop]][,paste0('original.max_',lock.loop)] <-
              max(abs(.original.mean[samples.to.compare]), na.rm = TRUE)
            elec.mean.maxes.pre.post[[patient]][[task.loop]][[elec.loop]][,paste0('original.peak.time_',lock.loop)] <- .original.peak.time
            
            rm(.original.mean)
        }; rm(lock.loop)
        } # if(.n.trials > 0)
      }; rm(elec.loop)
    }; rm(task.loop)

    ### Clean up
    elec.means.pre.post[[patient]] <- 
      bind_rows(lapply(elec.means.pre.post[[patient]], function(x){bind_rows(x)}))
    elec.mean.maxes.pre.post[[patient]] <-
      lapply(elec.mean.maxes.pre.post[[patient]], bind_rows)
    
    ### Remove the "keep" data from the "adjust" data
    rm(data, keep.data.stimulus, keep.data.production)
  }; rm(patient)
  
  
  ###
  ### Save data
  ###
  
  ### Save all data
  save.these <- c(
    'patient.trial.info',
    'median.rt.samples',
    'stretch.data',
    'stretch.sample.labels',
    'elec.mean.maxes.pre.post'
  )
  
  save.data.path <- paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',
                           freq.band.loop,
                           '/data/')
  dir.create(save.data.path, showWarnings = FALSE, recursive = TRUE)
  save(list = save.these,
       file = paste0(save.data.path,'elec data with RTs stretched to global median by task - bad trials and .025 to .95 RT outliers excluded.RData'))
  
  ### Save just median.rt.samples
  save.data.path <- paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/')
  dir.create(save.data.path, showWarnings = FALSE, recursive = TRUE)
  save(median.rt.samples,
       file = paste0(save.data.path,'median RT samples.RData'))
  
  ### Save just stretch sample labels
  save.data.path <- paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/')
  dir.create(save.data.path, showWarnings = FALSE, recursive = TRUE)
  save(stretch.sample.labels,
       file = paste0(save.data.path,'warped sample labels.RData'))
  
} # freq.band.loop

# Finish!
message('Script completed successfully. ',Sys.time())
