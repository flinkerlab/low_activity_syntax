### Perform task comparisons (e.g., sentence >? list) for each electrode at each time sample
### adam.milton.morgan@gmail.com

###
### Readme
###

###
### Setup
###

### Packages
library('beepr') # for playing beeps
library('pals')
library('viridis')
library('grid')
library('gridExtra')
library('doParallel') # for parallel processing
library('foreach') # for parallel for loops
library('BayesFactor') # for Bayes Factor t-tests
library('effsize') # for getting Cohen's d


# Clean up
rm(list=ls())
cat("\014")

# Patients
patients = c('Patient001',
             'Patient002',
             'Patient003', 
             'Patient004',
             'Patient005',
             # 'Patient006', # didn't do the List Production block 
             'Patient007',
             'Patient008', 
             'Patient009',
             'Patient010')

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 2
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 5
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 10
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 28
}

### Lemmas
# Function to convert between samples, sample labels, and times
source(paste0(path,'/analysis/R/functions/time_convert.R'))
# Function for rolling average smoothing
source(paste0(path,'/analysis/R/functions/smoothing.R'))
# Function for plotting stacked trials in z (color) dimension
source(paste0(path,'/analysis/R/functions/stack_plot.R'))
# Close any old parallel backends
source(paste0(path,'/analysis/R/functions/unregister_dopar.R'))
# Just good trial functions
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))


### Set up elec info
elec.info <- 
  read.csv(paste0(path,
                  'analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
elec.info$use.these.elecs <- 
  as.numeric((! elec.info$region_clinical %in% c('Unknown',
                                                 '',
                                                 'NaN',
                                                 'Left-Cerebral-White-Matter',
                                                 'Left-Inf-Lat-Vent',
                                                 'Right-Cerebral-White-Matter')) &
               (! is.na(elec.info$region_clinical)) &
               (! is.na(elec.info$MNI_x)) &
               (! is.na(elec.info$MNI_y)) &
               (! is.na(elec.info$MNI_z)) &
               (elec.info$bad_elec == 0)# &
             # (elec.info$visual_elec == 0) &
             # (elec.info$active == 1) & 
             # (!substr(elec.info$elec,1,1) == "D") # No depth elecs -- first letter of elec label "D"
  )
rownames(elec.info) <- elec.info$patient_elec
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec

### Loop thru frequency bands
keep <- c(ls(), 'keep', 'band.loop', 'break.at.median.loop')

# for(band.loop in c('high_gamma', 'beta')){ ## UNCOMMENT
band.loop = c('high_gamma', 'beta')[1]

# Two ways of splitting data - by median SP or median LP time
for(break.at.median.loop in c('sp','lp')[2]){ ## UNCOMMENT
  # break.at.median.loop = c('sp','lp')[2]
  
  # Clean up
  rm(list = ls()[! ls() %in% keep])
  gc()
  
  # Alphas
  alphas.to.try <- data.frame('alpha' = c(".05",".01",".001"),
                              'min_ms' = c(100, 100, 50))
  alphas.to.try$label <- paste0('alpha=',alphas.to.try$alpha,'_min=',alphas.to.try$min_ms,'ms')
  rownames(alphas.to.try) <- alphas.to.try$label
  
  message('Beginning stats on ',band.loop,' - break.at.median.loop="',break.at.median.loop,'". ',Sys.time())
  
  ### Loop thru patients in parallel
  ## Set up parallel processing
  unregister_dopar()
  n.cores <- n.cores.to.use
  cl <- makeCluster(n.cores, type="FORK")
  registerDoParallel(cl)
  
  ### Loop thru warped/unwarped datasets
  patient.stats <- foreach(patient = patients) %dopar% { ## UNCOMMENT
    # patient = patients[1]
    message("Beginning patient ",patient," data. ", Sys.time())
    
    
    ###
    ### Data and metadata
    ###
    
    # Define words
    words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')
    # Patient Patient008 took like twice as long to produce Drac and Frank and has artifactual differences between those words and others; remove
    if(patient == 'Patient008'){
      words <- c('chicken','dog','ninja','nurse')
    }
    
    ## Get data
    ### Load data
    message('Attaching electrode data...')
    attach(paste0(path,
                  'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',
                  band.loop,
                  '/data/elec data with RTs stretched to global median by task - bad trials and .025 to .95 RT outliers excluded.RData'))
    
    # Get data
    elec.data <- stretch.data[[patient]]
    trial.info <- patient.trial.info[[patient]]
    tasks <- names(elec.data)
    
    # Add blank dataframe for 837's listing block
    if(patient == 'Patient006'){
      blank.Patient006.lp.df <- data.frame(matrix(nrow = 0,
                                             ncol = length(stretch.sample.labels$lp[["locked_to_production_onset"]])))
      colnames(blank.Patient006.lp.df) <- stretch.sample.labels$lp[["locked_to_production_onset"]]
      for(elec.loop in names(elec.data$pn)){
        elec.data$lp[[elec.loop]] <- blank.Patient006.lp.df
      }; rm(elec.loop)
      rm(blank.Patient006.lp.df)
    }
    
    detach()
    message('...detached!')
    
    # Get median RTs
    load(paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
    median.rt.times <- time.convert(median.rt.samples, "samples", "times")
    
    # Sample labels
    sample.labels <- colnames(elec.data[['pn']][[1]])
    times <- time.convert(sample.labels, "sample.labels", "times")
    
    # Get electrode labels
    elec.labels <- names(elec.data[['pn']])
    
    # Subset to just active elecs with good localizations
    elec.labels <- elec.labels[elec.labels %in% use.these.elecs]
    
    
    ###
    ### Loop thru elecs and get mean HPG
    ###
    
    
    ### Get time limits for SP and LP:
    if(break.at.median.loop == 'sp'){
      # Align to stimulus and take the first median.rt['sp']/2
      # Align to production and take the last median.rt['sp']/2 until 400ms(= the latest time for which we've got RSIs)
      sp.time.boundaries.samples <- list('locked_to_stimulus_onset' = c('start' = -as.numeric(median.rt.samples['sp']), 'end' = floor(-(as.numeric(median.rt.samples['sp']) + .1)/2)),
                                         'locked_to_production_onset' = c('start' = ceiling(-(as.numeric(median.rt.samples['sp']) + .1)/2), 'end' = 400))
      half.median.n.samples <- as.numeric(sp.time.boundaries.samples[['locked_to_stimulus_onset']]['end']) - 
        as.numeric(sp.time.boundaries.samples[['locked_to_stimulus_onset']]['start'])
      lp.time.boundaries.samples <- list('locked_to_stimulus_onset' = c('start' = -as.numeric(median.rt.samples['lp']), 
                                                                        'end' = -as.numeric(median.rt.samples['lp']) + half.median.n.samples),
                                         'locked_to_production_onset' = sp.time.boundaries.samples[['locked_to_production_onset']])
      sp.time.boundaries <- lapply(sp.time.boundaries.samples, function(x){time.convert(x, "samples", "times")})
      lp.time.boundaries  <- lapply(lp.time.boundaries.samples, function(x){time.convert(x, "samples", "times")})  
    }
    if(break.at.median.loop == 'lp'){
      # Align to stimulus and take the first median.rt['lp']/2
      # Align to production and take the last median.rt['lp']/2 until 400ms(= the latest time for which we've got RSIs)
      lp.time.boundaries.samples <- list('locked_to_stimulus_onset' = c('start' = -as.numeric(median.rt.samples['lp']), 
                                                                        'end' = floor(-(as.numeric(median.rt.samples['lp']) + .1)/2)),
                                         'locked_to_production_onset' = c('start' = ceiling(-(as.numeric(median.rt.samples['lp']) + .1)/2), 
                                                                          'end' = 400))
      half.median.n.samples <- as.numeric(lp.time.boundaries.samples[['locked_to_stimulus_onset']]['end']) - 
        as.numeric(lp.time.boundaries.samples[['locked_to_stimulus_onset']]['start'])
      sp.time.boundaries.samples <- list('locked_to_stimulus_onset' = c('start' = -as.numeric(median.rt.samples['sp']), 
                                                                        'end' = -as.numeric(median.rt.samples['sp']) + half.median.n.samples),
                                         'locked_to_production_onset' = lp.time.boundaries.samples[['locked_to_production_onset']])
      
      # Pad by 1/2 min threshold to be fair here -- less possibility to be sig
      padding <- time.convert(50, "times", "samples")
      sp.time.boundaries.samples$locked_to_stimulus_onset['end'] <- sp.time.boundaries.samples$locked_to_stimulus_onset['end'] + padding
      sp.time.boundaries.samples$locked_to_production_onset['start'] <- sp.time.boundaries.samples$locked_to_production_onset['start'] - padding
      lp.time.boundaries.samples$locked_to_stimulus_onset['end'] <- lp.time.boundaries.samples$locked_to_stimulus_onset['end'] + padding
      lp.time.boundaries.samples$locked_to_production_onset['start'] <- lp.time.boundaries.samples$locked_to_production_onset['start'] - padding
      
      # Get boundaries in ms
      sp.time.boundaries <- lapply(sp.time.boundaries.samples, function(x){time.convert(x, "samples", "times")})
      lp.time.boundaries  <- lapply(lp.time.boundaries.samples, function(x){time.convert(x, "samples", "times")})
    }
    
    
    
    
    
    # Sample labels
    use.these.samples <- list('sp' = list(), 'lp' = list())
    for(lock.loop in names(sp.time.boundaries.samples)){
      use.these.samples[['sp']][[lock.loop]] <- 
        time.convert(sp.time.boundaries.samples[[lock.loop]][['start']]:sp.time.boundaries.samples[[lock.loop]][['end']], 
                     "samples", "sample.labels")
      use.these.samples[['lp']][[lock.loop]] <- 
        time.convert(lp.time.boundaries.samples[[lock.loop]][['start']]:lp.time.boundaries.samples[[lock.loop]][['end']], 
                     "samples", "sample.labels")
    }; rm(lock.loop)
    
    
    ### Loop thru samples and do comparisons
    sp.greater.than.lp.ps <- list()
    for(lock.loop in names(sp.time.boundaries.samples)){
      # lock.loop = names(sp.time.boundaries.samples)[1]
      
      sp.greater.than.lp.ps[[lock.loop]] <- data.frame(matrix(nrow = length(use.these.samples[['sp']][[lock.loop]]),
                                                              ncol = length(elec.labels),
                                                              dimnames = list(use.these.samples[['sp']][[lock.loop]],
                                                                              elec.labels)))
      
      for(sample.loop in 1:length(use.these.samples[['sp']][[lock.loop]])){
        # sample.loop = 1
        
        sp.current.sample <- use.these.samples[['sp']][[lock.loop]][sample.loop]
        lp.current.sample <- use.these.samples[['lp']][[lock.loop]][sample.loop]
        
        sp.current.data <- lapply(elec.data[['sp']][elec.labels], function(x){x[,sp.current.sample]})
        lp.current.data <- lapply(elec.data[['lp']][elec.labels], function(x){x[,lp.current.sample]})
        
        for(elec.loop in names(sp.current.data)){
          # elec.loop = names(sp.current.data)[1]
          
          sp.greater.than.lp.ps[[lock.loop]][sp.current.sample, elec.loop] <- 
            wilcox.test(sp.current.data[[elec.loop]],
                        lp.current.data[[elec.loop]], 
                        alternative = "greater")$p.value
          
        }; rm(elec.loop)
        # message('Sample loop ',sample.loop,' of ',length(use.these.samples[['sp']][[lock.loop]]),' complete. ',Sys.time())
      }; rm(sample.loop)
    }; rm(lock.loop)
    
    
    ### Get sig
    sp.greater.than.lp.sig <- list()
    sp.greater.than.lp.sig.windows <- list()
    sp.greater.than.lp.sig.elecs <- list()
    for(alpha.loop in 1:nrow(alphas.to.try)){
      # alpha.loop = alphas.to.try[1]
      
      current.alpha = as.numeric(alphas.to.try$alpha[alpha.loop])
      current.min.ms = alphas.to.try$min_ms[alpha.loop]
      alpha.label <- alphas.to.try$label[alpha.loop]
      
      sp.greater.than.lp.sig[[alpha.label]] <- list()
      sp.greater.than.lp.sig.windows[[alpha.label]] <- list()
      for(lock.loop in names(sp.greater.than.lp.ps)){
        # lock.loop = names(sp.greater.than.lp.ps)[1]
        
        sp.greater.than.lp.sig[[alpha.label]][[lock.loop]] <- apply(sp.greater.than.lp.ps[[lock.loop]], 2, function(x){as.numeric(x < current.alpha)})
        rownames(sp.greater.than.lp.sig[[alpha.label]][[lock.loop]]) <- rownames(sp.greater.than.lp.ps[[lock.loop]])
        
        sp.greater.than.lp.sig.windows[[alpha.label]][[lock.loop]] <- apply(sp.greater.than.lp.sig[[alpha.label]][[lock.loop]], 2, function(x){
          x <- get.significant.windows(x,
                                       .sample.labels = rownames(sp.greater.than.lp.ps[[lock.loop]]),
                                       output.class = "data.frame",
                                       include.duration = TRUE,
                                       .exclude.sig.durations.under.ms = current.min.ms)
          return(x)
        })
        
      }; rm(lock.loop)
      
      # Combine sig windows
      for(elec.loop in elec.labels){
        # elec.loop = elec.labels[1]
        sp.greater.than.lp.sig.windows[[alpha.label]][[elec.loop]] <-
          rbind(sp.greater.than.lp.sig.windows[[alpha.label]]$locked_to_stimulus_onset[[elec.loop]],
                sp.greater.than.lp.sig.windows[[alpha.label]]$locked_to_production_onset[[elec.loop]])
      }; rm(elec.loop)
      
      # Clean up
      sp.greater.than.lp.sig.windows[[alpha.label]]$locked_to_stimulus_onset <- NULL
      sp.greater.than.lp.sig.windows[[alpha.label]]$locked_to_production_onset <- NULL
      
      # Get sig elecs
      sp.greater.than.lp.sig.elecs[[alpha.label]] <- names(which(sapply(sp.greater.than.lp.sig.windows[[alpha.label]], function(x){nrow(x) > 0})))
      
    }; rm(alpha.loop, alpha.label)
    
    
    ###
    ### Save
    ###
    
    # current.save.dir <- paste0(path,
    #                            'analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/data/',
    #                            'warped data/',
    #                            band.loop,
    #                            '/patients - individually/',
    #                            'from stimulus to 400ms post speech onset/')
    # dir.create(current.save.dir, showWarnings = FALSE, recursive = TRUE)
    # save.these <- c('sp.greater.than.lp.ps',
    #                 'sp.greater.than.lp.sig.windows',
    #                 'sp.greater.than.lp.sig.elecs')
    # save(list = save.these,
    #      file = paste0(current.save.dir,patient,' - sp vs. lp comparisons - all time.RData'))
    
    # Don't return anything from foreach loop so as not to bog down memory
    output <- list('sp.greater.than.lp.ps' = sp.greater.than.lp.ps,
                   'sp.greater.than.lp.sig.windows' = sp.greater.than.lp.sig.windows,
                   'sp.greater.than.lp.sig.elecs' = sp.greater.than.lp.sig.elecs)
    return(output)
    
  } # patient, lock.loop, warp.loop
  # End parallel processing
  stopCluster(cl)
  
  beep()  
  
  
  ### Reorganize
  sp.greater.than.lp.ps <- lapply(patient.stats, function(x){x[['sp.greater.than.lp.ps']]})
  sp.greater.than.lp.sig.windows.temp <- lapply(patient.stats, function(x){x[['sp.greater.than.lp.sig.windows']]})
  sp.greater.than.lp.sig.elecs.temp <- lapply(patient.stats, function(x){x[['sp.greater.than.lp.sig.elecs']]})
  
  
  
  ### Collapse across patients
  sp.greater.lp.sig.windows <- list()
  sp.greater.lp.sig.elecs <- list()
  for(alpha.loop in alphas.to.try$label){
    # alpha.loop = alphas.to.try$label[1]
    
    # Combine sig elecs
    sp.greater.lp.sig.elecs[[alpha.loop]] <- unlist(lapply(sp.greater.than.lp.sig.elecs.temp, function(x){x[[alpha.loop]]}))
    
    # Combine sig windows
    sp.greater.lp.sig.windows[[alpha.loop]] <- list()
    for(patient.loop in 1:length(sp.greater.than.lp.sig.windows.temp)){
      # patient.loop = 1
      
      for(elec.loop in names(sp.greater.than.lp.sig.windows.temp[[patient.loop]][[alpha.loop]])){
        # elec.loop = names(sp.greater.than.lp.sig.windows.temp[[patient.loop]][[alpha.loop]])[1]
        
        sp.greater.lp.sig.windows[[alpha.loop]][[elec.loop]] <-
          sp.greater.than.lp.sig.windows.temp[[patient.loop]][[alpha.loop]][[elec.loop]]
        
      }; rm(elec.loop)
      
    }; rm(patient.loop)
    
  }; rm(alpha.loop)
  
  
  ### Save
  
  current.save.dir <- paste0(path,
                             'analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/data/',
                             'warped data/',
                             band.loop,
                             '/patients - combined/',
                             'from stimulus to 400ms post speech onset/')
  dir.create(current.save.dir, showWarnings = FALSE, recursive = TRUE)
  save.these <- c('sp.greater.than.lp.ps',
                  'sp.greater.lp.sig.elecs',
                  'sp.greater.lp.sig.windows')
  save(list = save.these,
       file = paste0(current.save.dir,'sp vs. lp comparisons - all times to and from median rt of ',break.at.median.loop,'.RData'))
  
}#; rm(break.at.median.loop)









# } # band.loop


















message('Script complete. ', Sys.time())

