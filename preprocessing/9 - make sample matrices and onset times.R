### Chicken Syntax
### Pre-processing stage 9
### Make sample epochs and remove bad trials
### July 2020

## This whole script is redundant and at some point should be made to be part of preprocessing stage 7, but going for quick and easy right now

rm(list=ls())

### Set paths
if(Sys.info()['sysname'] == 'Darwin'){
  path = '/Users/am4611/Dropbox/Research/ChickenSyntax/data/'
  n.cores.to.use <- 8
} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'
  n.cores.to.use <- 24
} # Ubuntu

#####################
### Manually edit ###
#####################
patient <- 'Patient010'
srate <- 512 # sampling rate
time.lock.options <- c('locked_to_production_onset','locked_to_stimulus_onset')
smoothing.options <- 'multiband_analytic_amplitude' #c('unsmoothed','smoothed')
#n.smoothing.samples.post.and.pre = 40
type.of.baselining <- c('z_scored','z_scored_beta','raw_potentials_postCAR')
(start.start.time <- Sys.time())

###################
### Auto - rest ###
###################
## Packages
library("tidyverse")
library("beepr")
library('doParallel') # for parallel processing
library('foreach') # for parallel for loops
# Close any old parallel backends
source(paste0(path,'../analysis/R/functions/unregister_dopar.R'))


# Where's the raw data?
epoched.data.dir <- paste0(path,patient,'/data/epoched/')


### Metadata
# Trials to remove (epileptic activity)
bad.trials <- read.csv(paste0(epoched.data.dir,patient,'_trials_to_remove.csv'))$trials_to_remove

# Read in number of trials prior to trial exclusions
nTrials_before_exclusions <- read.csv(paste0(epoched.data.dir,patient,'_nTrials_before_trial_exclusions.csv'),header=FALSE,col.names='n')$n

# Read in trial labels (dog, chicken, ..., eventually also subject/object, DP/NP)
trial.labels <- read.csv(paste0(epoched.data.dir,patient,
                                '_trial_labels_table.csv'))
# Replace blanks with NAs
for(col in 1:ncol(trial.labels)){
  trial.labels[which(trial.labels[,col]==""),col] <- NA
  if(! is.numeric(trial.labels[,col])){trial.labels[,col] <- factor(trial.labels[,col])}
}; rm(col)

# Remove bad trials
trial.labels <- trial.labels[-bad.trials,]
# Write trial labels without bad trials to file
write.csv(trial.labels,
          paste0(epoched.data.dir,patient,'_trial_labels_without_bad_trials.csv'),
          row.names=FALSE, quote=FALSE)

# Number of trial label columns (to be removed for averaging, etc)
n.lx.cols <- ncol(trial.labels)
write.table(n.lx.cols,
            paste0(epoched.data.dir,patient,'_n_linguistic_columns.txt'),
            row.names=FALSE, col.names=FALSE, quote=FALSE)


## Set up parallel processing
# Close any old parallel backends
unregister_dopar()
# Set up parallel workers in case caret wants to use it
cl <- makeCluster(n.cores.to.use, type="FORK")
registerDoParallel(cl)

### Loop through smoothing and time-lock options:
foreach(smooth.loop = smoothing.options) %:%
  # (smooth.loop = smoothing.options[1])
  
  ## Now loop thru time lock options
foreach(time.lock = time.lock.options) %:%
    # (time.lock = time.lock.options[1])
    
  foreach(baseline.loop = type.of.baselining) %dopar% {
      # (baseline.loop = type.of.baselining[1])
      message(patient,": Starting '",smooth.loop,"', '",time.lock,"', and '",baseline.loop,"' loops. ", Sys.time())
      
      ### Load metadata
      # Directories
      elec.data.dir <- paste0(epoched.data.dir,
                              '/electrode_data/',
                              smooth.loop,'/',
                              time.lock,'/',
                              baseline.loop,'/')
      
      # What are the names of the files (i.e., names of electrodes)?
      electrodes <- gsub(".csv","",list.files(elec.data.dir))
      
      ## Read in sample labels (e.g., "sample_405_post")
      if(time.lock == "locked_to_production_onset"){
        sample.labels <- as.character(read.csv(paste0(epoched.data.dir,patient,'_sample_labels.csv'),
                                               header=FALSE,
                                               col.names='sample')$sample)
        # Get earliest sample (absolute value)
        earliest.time <- 1201 # the first time you'll want to look at pre-onset; always positive; should be less than the absolute earliest number of samples available in the electrode data plus the window.duration that you're going to average any values over. so if you've got electrode epochs starting at 1200ms and you want to look from 1000ms but you're averaging over 250ms window's that's gonna be a problem (it's going to try to average the HGP from 1250ms to 1000ms pre onset and call that 1000ms)
        earliest.sample <- round(earliest.time * srate / 1000, 0) # as.numeric(gsub("sample_","",sample.labels[1]))
        latest.time <- 701 # the last time you'll want to look at, must be <= greatest post-onset time available; always positive
        latest.sample <- round(latest.time * srate / 1000, 0)
        rm(earliest.time, latest.time)
        
        # Translate sample number into time and label 
        samples.df <- data.frame('sample'=(-earliest.sample):latest.sample)
        samples.df$time <- round(samples.df$sample * 1000 / srate, 0) 
        samples.df$label <- gsub("-","neg",paste0("sample_",samples.df$sample))
      }
      if(time.lock == "locked_to_stimulus_onset"){
        sample.labels <- as.character(read.csv(paste0(epoched.data.dir,patient,'_sample_labels_stimLocked.csv'),
                                               header=FALSE,
                                               col.names='sample')$sample)
        earliest.time <- 701 
        earliest.sample <- round(earliest.time * srate / 1000, 0) 
        latest.time <- 1201 
        latest.sample <- round(latest.time * srate / 1000, 0)
        rm(earliest.time, latest.time)
        
        # Translate sample number into time and label 
        samples.df <- data.frame('sample'=(-earliest.sample):latest.sample)
        samples.df$time <- round(samples.df$sample * 1000 / srate, 0) 
        samples.df$label <- gsub("-","neg",paste0("sample_",samples.df$sample))
      }
      
      ## Read in all data
      all.data <- list()
      for(electrode.loop in 1:length(electrodes)){ # electrode.loop = 1 # for troubleshooting
        current.electrode <- electrodes[electrode.loop]
        # Read data
        elec.data <- read.csv(paste0(elec.data.dir,current.electrode,'.csv'),
                              header=FALSE,
                              col.names=sample.labels)
        
        # If bad trials *have* already been removed, then the top row is probably a duplicate of the col.names. Remove it.
        if(!elec.data[1,1] == gsub("sample_","",elec.data[1,1])){ # all(elec.data[1,] == colnames(elec.data))
          elec.data <- elec.data[-1,]
        }
        # If bad trials haven't already been removed (as in, this script has been run already), remove them
        if(nrow(elec.data) == nTrials_before_exclusions){
          elec.data <- elec.data[-bad.trials,]
        }
        
        all.data[[current.electrode]] <- cbind(trial.labels,
                                               elec.data)
        # Progress update and clean up
        if(electrode.loop %% 50 == 0){
          message("Reading in elec ", electrode.loop, " of ",length(electrodes)," complete. ",Sys.time())}
        rm(elec.data, current.electrode)
      } # takes ~5 mins
      # beep(2);Sys.sleep(3);beep(3);Sys.sleep(3);beep(4);Sys.sleep(3);beep(5);Sys.sleep(3);beep(6)
      
      ## Store the electrode data (with bad trials removed) as .Rdata file (compressed)
      elec.data.save.path <- paste0(path,
                                    patient,
                                    '/data/epoched/elec_data_without_bad_trials/',
                                    smooth.loop,'/',
                                    time.lock,'/',
                                    baseline.loop,'/')
      dir.create(elec.data.save.path, showWarnings = FALSE, recursive = TRUE)
      save(list = c('all.data'),
           file = paste0(elec.data.save.path,
                         patient,
                         ' elec data.RData'))
      
      
      ## Loop through samples, creating a list of dataframes for each window
      # Each dataframe should be nTrial rows by nElectrode columns
      sample.data <- list()
      for(sample.loop in 1:nrow(samples.df)){ # sample.loop = 1
        current.sample <- samples.df$sample[sample.loop]
        current.sample.label <- samples.df$label[sample.loop]
        
        # Initialize dataframe
        sample.data[[current.sample.label]] <- trial.labels
        for(electrode.loop in 1:length(electrodes)){ # electrode.loop <- 1 # for troubleshooting
          current.electrode <- electrodes[electrode.loop]
          sample.data[[current.sample.label]][,current.electrode] <- all.data[[current.electrode]][,current.sample.label]
        }
        
        if(sample.loop %% 100 == 0){
          message("Creating sample ", sample.loop, " of ",nrow(samples.df),". ",Sys.time())}
      } # sample.loop
      
      # Create save directory if doesn't exist
      dir.create(paste0(epoched.data.dir,
                        '/sample_data_512Hz_without_bad_trials/',
                        smooth.loop,'/',
                        time.lock,'/',
                        baseline.loop,'/'),
                 showWarnings = FALSE,
                 recursive = TRUE)
      
      # Save
      save(list = c('sample.data'),
           file = paste0(epoched.data.dir,
                         '/sample_data_512Hz_without_bad_trials/',
                         smooth.loop,'/',
                         time.lock,'/',
                         baseline.loop,'/',
                         patient,
                         ' sample data.RData'))
      
      # Write window labels to file
      write.csv(samples.df,
                paste0(epoched.data.dir,
                       '/sample_data_512Hz_without_bad_trials/',
                       smooth.loop,'/',
                       time.lock,'/',
                       patient,
                       '_sample_time_translation.csv'),
                row.names=FALSE, quote=FALSE)
    # } # baseline.loop
  # } # time.lock
} # smooth.loop

# End parallel processing
stopCluster(cl)
unregister_dopar()
beep()

# Delete the original electrode data file (with all of the bad trials)
# unlink(paste0(epoched.data.dir,'electrode_data/'),
#        recursive = TRUE)
end.end.time <- Sys.time()
message("Total script duration: ",round(end.end.time - start.start.time, 1)) # 25 minutes once, 2.6 hours another time, both on Lambda...

