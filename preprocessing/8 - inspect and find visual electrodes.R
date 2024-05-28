### Chicken Syntax
### Pre-processing stage 8
### Make window epochs
### July 2020

## This whole script is redundant and at some point should be made to be part of preprocessing stage 7, but going for quick and easy right now

rm(list=ls())

####################
### Manully edit ###
####################
patient <- 'Patient010'
srate <- 512 # sampling rate
type.of.baselining <- 'z_scored' # Folder names: z_scored, percent_change_baseline, hgp

#############################################
### Auto - more manual edits needed below ###
#############################################

## Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/am4611/Dropbox/Research/ChickenSyntax/data/'
  n.cores.to.use <- 4
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'
  n.cores.to.use <- 48
}

## Packages
library("tidyverse")
library('beepr')
library('doParallel') # for parallel processing
library('foreach') # for parallel for loops
# Close any old parallel backends
source(paste0(path,'../analysis/R/functions/unregister_dopar.R'))
source(paste0(path,'../analysis/R/functions/time_convert.R'))

### Load metadata
# Where are the files stored?
data.dir <- paste0(path,patient,'/data/epoched/electrode_data')
# What are the names of the files (i.e., names of electrodes)?
electrodes <- gsub(".csv","",list.files(paste0(data.dir,
                                               '/multiband_analytic_amplitude/locked_to_production_onset/',
                                               type.of.baselining)))
# Read in trial labels (dog, chicken, ..., eventually also subject/object, DP/NP)
trial.labels <- read.csv(paste0(data.dir,'/../',patient,'_trial_labels_table.csv'))
# Replace blanks with NAs
for(col in 1:ncol(trial.labels)){
  trial.labels[which(trial.labels[,col]==""),col] <- NA
  if(! is.numeric(trial.labels[,col])){trial.labels[,col] <- factor(trial.labels[,col])}
}; rm(col)

# Number of trial label columns (to be removed for averaging, etc)
n.lx.cols <- ncol(trial.labels)
# Which rows of data to use for classification (naming trials only)
naming.trials <- which(trial.labels$case == "none")
# What are the unique words? (dog, chicken, ...)
words <- as.character(unique(droplevels(subset(trial.labels, case=='none'))$word))


## Read in sample labels (e.g., "sample_405_post")
sample.labels <- list()
sample.labels[['locked_to_production_onset']] <-
  as.character(read.csv(paste0(data.dir,'/../',
                               patient,'_sample_labels.csv'),
                        header=FALSE,
                        col.names='sample')$sample)
sample.labels[['locked_to_stimulus_onset']] <- 
  as.character(read.csv(paste0(data.dir,'/../',
                               patient,'_sample_labels_stimLocked.csv'),
                        header=FALSE,
                        col.names='sample')$sample)

## Read in ROIs
localizations.df <- read.csv(paste(data.dir,'/../../../localization/',patient,'_electrode_localization.csv',sep=''))[,c('labels','roi_electrodes','T1_AnatomicalRegion')]
rois <- as.character(unique(localizations.df$roi_electrodes))
localizations.df$electrode_labels <- paste0(patient,
                                            "_",
                                            localizations.df$labels)
roi.electrodes <- list()
for(roi.loop in rois){ # roi.loop <- rois[3]
  roi.electrodes[[roi.loop]] <- as.character(subset(localizations.df, roi_electrodes == roi.loop)$electrode_labels)
}

### Read in all data

## Set up parallel processing
# Close any old parallel backends
unregister_dopar()
# Set up parallel workers in case caret wants to use it
cl <- makeCluster(n.cores.to.use, type="FORK")
registerDoParallel(cl)

## Read in all data
all.dat <- list()
for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
  all.dat[[lock.loop]] <- 
    foreach(electrode.loop = electrodes, 
            .options.nws = list(chunkSize = 5),
            .final = function(x) setNames(x, electrodes)) %dopar%{
      
      # Read data (and pre-onset cols) 
      elec.dat <- read.csv(paste0(data.dir,'/',
                                  'multiband_analytic_amplitude/',
                                  lock.loop,'/',
                                  type.of.baselining,'/',
                                  electrode.loop,
                                  '.csv'),
                           header=FALSE,
                           col.names=sample.labels[[lock.loop]])
      
      # Add linguistic columns
      elec.dat <- cbind(trial.labels, elec.dat)
    
      return(elec.dat)
  } # foreach
  message(lock.loop," loop completed at ",Sys.time())
} # ~15 mins

# End parallel processing
stopCluster(cl)
beep()


# Get row numbers for naming trials, and quick naming
naming.trial.rows <- with(trial.labels, 
                          which(case=='none' & is.na(quick_image_duration)))
naming.200.rows <- with(trial.labels,
                        which(quick_image_duration==200))
naming.300.rows <- with(trial.labels,
                        which(quick_image_duration==300))

## Define smoothing function
smoothing <- function(time.series, 
                      n.samples.pre = 20, 
                      n.samples.post = 20,
                      na.pad = TRUE){
  out <- c()
  for(i in (n.samples.pre + 1):(length(time.series)-n.samples.post)){
    out <- c(out,
             mean(time.series[(i-n.samples.pre):(i+n.samples.post)]))}
  if(na.pad){
    out <- c(rep(NA, times=n.samples.pre),
             out,
             rep(NA, times=n.samples.post))}
  return(out)
}

## Define function for plotting
mean.electrode.plot <- function(electrodes=electrodes,
                                lock.loop="locked_to_stimulus_onset"){
  for(electrode.loop in 1:length(electrodes)){ 
    # electrode.loop <- 1 # for troubleshooting
    current.electrode <- electrodes[electrode.loop]
    current.localization <- as.character(localizations.df[localizations.df$electrode_labels==current.electrode,'T1_AnatomicalRegion'])
    current.data <- all.dat[[lock.loop]][[current.electrode]]
    
    x <- as.numeric(gsub("sample_","",gsub("neg","-",sample.labels[[lock.loop]])))/512*1000
    y <- colMeans(current.data[naming.trial.rows,
                               (n.lx.cols+1):ncol(current.data)])
    y2 <- colMeans(current.data[naming.200.rows,
                                (n.lx.cols+1):ncol(current.data)])
    y3 <- colMeans(current.data[naming.300.rows,
                                (n.lx.cols+1):ncol(current.data)])
    plot(x=x,
         y=smoothing(y),
         type='l',
         lwd=20,
         main=paste0(patient," - ",
                     gsub(paste0(patient,"_"),"",current.electrode)," - ",
                     current.localization),
         ylab=type.of.baselining,
         xlab='time (ms)',
         col='lightgrey',
         ylim=c(min(smoothing(c(y,y2,y3), na.pad=FALSE)),
                max(smoothing(c(y,y2,y3), na.pad=FALSE))))
    lines(x=x,
          y=smoothing(y2),
          lwd=3, col='royalblue')
    lines(x=x,
          y=smoothing(y3),
          lwd=3, col='darkred')
    abline(v=0)
    if(lock.loop=='locked_to_stimulus_onset'){
      abline(v=200, lty=2, col='royalblue')
      abline(v=300, lty=2, col='darkred')}
    readline(prompt="Press [enter] to continue")
  }
}


## Make list of visual electrodes
# These should peak within 100ms or so of stimulus onset
# Compare with Leyao's functional mapping results to make sure they don't look radically different
# Look for visual electrodes and bad electrodes and add them below
for(lock.loop in c('locked_to_stimulus_onset','locked_to_production_onset')){
  # lock.loop == 'locked_to_stimulus_onset' # for troubleshooting
  message("Now beginning ",lock.loop,' plots.')
  readline(prompt="Press [enter] twice to continue")
  readline(prompt="(Nokh a mol!)")
  
  # Plot
  mean.electrode.plot(electrodes = electrodes,
                      lock.loop = lock.loop)
}

#####################
### Manually edit ###
#####################
# Add visual electrodes here (format example: 'G35')
visual.elecs <- c('G12','G18','G4') # Patient010 # maybe also 'G10','G2'
  
  

#############################################
### Auto - more manual edits needed below ###
#############################################
# Review visual and bad elecs
mean.electrode.plot(electrodes = paste0(patient,'_',visual.elecs))
# See these electrodes production-locked
mean.electrode.plot(electrodes = paste0(patient,'_',visual.elecs), 
                    lock.loop = 'locked_to_production_onset')

### Find outlier trials
## Maybe use raw potentials (all.dat.car) instead of HGP
## Function for plotting individual outlier trials
# Outlier defined as having the max value be >= the 95th %ile of trials' max values
bad.trials <- list()
plot.on <- FALSE
for(electrode.loop in 1:length(electrodes)){
  current.electrode <- electrodes[electrode.loop]
  current.localization <- as.character(localizations.df[localizations.df$electrode_labels==current.electrode,'T1_AnatomicalRegion'])
  current.data <- all.dat[['locked_to_production_onset']][[current.electrode]]
  
  # Find rows with big deviations
  maxes <- apply(current.data[,(n.lx.cols+1):ncol(current.data)],1,max)
  bad.trials[[current.electrode]] <- which(maxes > quantile(maxes,.95))
  
  # Plot
  x <- as.numeric(gsub("sample_","",gsub("neg","-",sample.labels[[lock.loop]])))/512*1000
  if(plot.on){
    plot(x=x,
         y=rep(NA,times=length(x)),
         type='l',
         ylim=c(-1,length(bad.trials[[current.electrode]])+2),
         main=paste0(patient," - ",
                     gsub(paste0(patient,"_"),"",current.electrode)," - ",
                     current.localization),
         ylab=type.of.baselining,
         xlab='time (ms)')
    abline(v=0, lty=2, col='grey')
    scaling.factor <- median(apply(current.data[bad.trials[[current.electrode]],(n.lx.cols+1):ncol(current.data)],1,max))
    for(i in 1:length(bad.trials[[current.electrode]])){
      lines(x=x,
            y=i + (current.data[bad.trials[[current.electrode]][i],(n.lx.cols+1):ncol(current.data)])/scaling.factor)
      text(x=x[1],
           y=i,
           labels=bad.trials[[current.electrode]][i],
           col='orange')
    }
    readline(prompt="Press [enter] to continue")
  }
}

# Remove visual elecs from list of potential bad trials
for(electrode.loop in visual.elecs){
  bad.trials[[electrode.loop]] <- NULL
}

# Find trials that crop up over and over again
bad.trials.t <- sort(table(factor(unlist(bad.trials))))
hist(bad.trials.t, breaks=40)
bad.trials.t.top <- bad.trials.t[which(bad.trials.t > quantile(bad.trials.t,.95))]
potential.exclude.trials <- sort(as.numeric(labels(bad.trials.t.top)[[1]]))

# Get the standard deviation for each electrode to use as scaling factor
electrode.sds <- list()
for(electrode.loop in 1:length(electrodes)){ # electrode.loop=1
  current.electrode <- electrodes[electrode.loop]
  current.data <- all.dat[['locked_to_production_onset']][[current.electrode]][,-(1:n.lx.cols)]
  electrode.sds[[current.electrode]] <- sd(unlist(current.data))
}

## Plot the bad trials for a few random electrodes
# Save plots
save.dir <- paste0(path,patient,'/preprocessing/stage 8 output - multiband/')
dir.create(save.dir)
par(mfrow = c(1,1))
n.electrodes.to.sample <- 45
lock.loop <- 'locked_to_production_onset'
x <- as.numeric(gsub("sample_","",gsub("neg","-",sample.labels[[lock.loop]])))/512*1000
plot.limits <- 512:1536
x <- x[plot.limits]
for(trial.loop in 1:length(potential.exclude.trials)){
  jpeg(paste0(save.dir,
              patient,
              "_Trial ",
              potential.exclude.trials[trial.loop],
              '.jpeg'),
       width = 512,
       height = 512)
  plot(x=x,
       y=rep(NA,times=length(x)),
       xlim=c(x[1]-500,x[length(x)]),
       ylim=c(-1,n.electrodes.to.sample+2),
       xlab=paste0("time (ms; ",lock.loop,")"),
       ylab="",
       main=paste0('Patient ',
                   patient,
                   ", Trial #",
                   potential.exclude.trials[trial.loop],
                   " (in top 95th %-ile ",
                   as.character(bad.trials.t.top[as.character(potential.exclude.trials[trial.loop])]),
                   ' times): ',
                   all.dat[[lock.loop]][[1]]$pos[potential.exclude.trials[trial.loop]],'-',
                   all.dat[[lock.loop]][[1]]$case[potential.exclude.trials[trial.loop]]))
  abline(v=c(-1000,0), lty=2, col='coral')
  abline(v= time.convert(-all.dat[[lock.loop]][[1]]$production_latency[potential.exclude.trials[trial.loop]], "samples", "times"), lty=2, col='forestgreen')
  sampled.electrodes <- sample(electrodes[which(!electrodes %in% paste0(patient,"_",visual.elecs))], n.electrodes.to.sample)
  for(electrode.loop in 1:n.electrodes.to.sample){
    scaling.factor <- electrode.sds[[sampled.electrodes[electrode.loop]]] * 10
    current.data <- all.dat[[lock.loop]][[sampled.electrodes[electrode.loop]]][potential.exclude.trials[trial.loop],]
    current.data <- current.data[,(n.lx.cols+1):ncol(current.data)]
    lines(x=x,
          y=(electrode.loop + current.data/scaling.factor)[plot.limits])
    text(x=x[1]-300,
         y=electrode.loop,
         labels=sampled.electrodes[electrode.loop],
         col='coral')
  }
  dev.off()
  message('Patient ',patient,
          ", Trial #",potential.exclude.trials[trial.loop],
          "(",trial.loop,
          " of ",length(potential.exclude.trials),") complete. ",
          Sys.time())
}
beep()
  
#####################
### Manually edit ###
#####################
# Store trials with epileptic activity
remove.trials <- c(34,123,144,188,319,440,441,708,709,772,870,871,878,879,976,1038,1044,1058,1074) # Patient010 new as of 7/10/2023


###################
### Auto - rest ###
###################

# Add back in any wrong_noun trials
susbtitution.trials <- which(! is.na(all.dat[[1]][[1]]$substitution_error))
add.trials.back <- susbtitution.trials[which(susbtitution.trials %in% remove.trials)]
remove.trials <- remove.trials[which(! remove.trials %in% add.trials.back)]


# Write these to file
write.csv(data.frame('visual_elecs'=visual.elecs),
          paste0(data.dir,'/../../',patient,'_visual_elecs.csv'),
          row.names=FALSE,quote=FALSE)
write.csv(data.frame('trials_to_remove'=remove.trials),
          paste0(data.dir,'/../',patient,'_trials_to_remove.csv'),
          row.names=FALSE,quote=FALSE)

