### Perform task comparisons (e.g., sentence >? list) for each electrode, in time bins (e.g., 250ms windows)
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
             'Patient006', 
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
    n.cores.to.use = 6
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 12
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
               (elec.info$bad_elec == 0) &
               # (elec.info$visual_elec == 0) &
               (elec.info$active == 1) & 
               (!substr(elec.info$elec,1,1) == "D") # No depth elecs -- first letter of elec label "D"
             )
rownames(elec.info) <- elec.info$patient_elec
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec

# Alphas
alphas.to.try <- c(.05,.01,.001)

### Loop thru frequency bands
keep <- c(ls(), 'keep', 'band.loop')
for(band.loop in c('high_gamma', 'beta')){
  # band.loop = c('high_gamma', 'beta')[1]
  rm(list = ls()[! ls() %in% keep])
  
  message('Beginning stats on ',band.loop,'. ',Sys.time())
  
  ### Loop thru patients in parallel
  ## Set up parallel processing
  unregister_dopar()
  n.cores <- n.cores.to.use
  cl <- makeCluster(n.cores, type="FORK")
  registerDoParallel(cl)
  
  ### Loop thru warped/unwarped datasets
  foreach(warp.loop = c('warped','unwarped')) %:%
    # warp.loop = c('warped','unwarped')[1]
    
    foreach(lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')) %:%
    # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[2]
    
    foreach(patient = patients) %dopar% {
      # patient = patients[1]
      message("Beginning patient ",patient,", ",lock.loop,", ",warp.loop," data. ", Sys.time())
      
      
      ###
      ### Data and metadata
      ###
      
      # Define words
      words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')
      # Patient008 took like twice as long to produce Drac and Frank and has artifactual differences between those words and others; remove
      if(patient == 'Patient008'){
        words <- c('chicken','dog','ninja','nurse')
      }
      
      ## Read in electrode data
      if(warp.loop == 'warped'){
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
                                                 ncol = length(stretch.sample.labels$lp[[lock.loop]])))
          colnames(blank.Patient006.lp.df) <- stretch.sample.labels$lp[[lock.loop]]
          for(elec.loop in names(elec.data$pn)){
            elec.data$lp[[elec.loop]] <- blank.Patient006.lp.df
          }; rm(elec.loop)
          rm(blank.Patient006.lp.df)
        }
        
        # Change column labels if locked to stimulus onset
        if(lock.loop == 'locked_to_stimulus_onset'){
          for(task.loop in tasks){
            # task.loop = tasks[1]
            elec.data[[task.loop]] <- lapply(elec.data[[task.loop]], function(x){
              colnames(x) <- stretch.sample.labels[[task.loop]][[lock.loop]]
              return(x)
            })
          }
        } # if lock.loop
        detach()
        message('...detached!')
        
      }else{ # if(warp.loop == 'warped'){}
        
        ## Get data
        elec.data <- list()
        message('Attaching electrode data...')
        attach(paste0(path,
                      'data/',
                      patient,
                      '/data/epoched/elec_data_without_bad_trials/multiband/',
                      lock.loop,
                      ifelse(band.loop == 'high_gamma','/z_scored/','/z_scored_beta/'),
                      patient,' elec data.RData'))
        
        # ECoG data
        elec.data[['pn']] <- lapply(all.data, just.good.pic.naming, .words = words)
        elec.data[['sp']] <- lapply(all.data, just.first.word.sentence.data)
        elec.data[['lp']] <- lapply(all.data, just.first.word.list.data)
        tasks <- names(elec.data)
        detach()
        message('...detached!')
        
        ## Number of linguistic features (columns)
        n.lx.cols = read.table(paste0(path,'data/',
                                      patient, 
                                      '/data/epoched/',
                                      patient,'_n_linguistic_columns.txt'))[1,1]
        
        # Remove linguistic columns
        for(task.loop in tasks){
          elec.data[[task.loop]] <- lapply(elec.data[[task.loop]], function(x){
            x <- x[,-c(1:n.lx.cols)]
            return(x)
          })
        }; rm(task.loop)
        
        ## Read in trial labels (linguistic data)
        trial.info.all <- read.csv(paste0(path,'data/',patient,'/data/epoched/',patient,'_trial_labels_without_bad_trials.csv'))
        # Factorize
        for(lx.col.loop in 1:ncol(trial.info.all)){ # lx.col.loop = 1
          if(class(trial.info.all[,lx.col.loop])=="character"){
            trial.info.all[,lx.col.loop] <- factor(trial.info.all[,lx.col.loop])}
        }; rm(lx.col.loop)
        
        # Split
        trial.info <- list()
        trial.info[['pn']] <- just.good.pic.naming(trial.info.all, .words = words)
        trial.info[['sp']] <- just.first.word.sentence.data(trial.info.all)
        trial.info[['lp']] <- just.first.word.list.data(trial.info.all)
        rm(trial.info.all)
        
      } # if(warp.loop == 'warped'){}else{
      gc()
      
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
      
      ## Average across bins of various sizes
      for(bin.loop in c(250,100,150,200,50)){ # in milliseconds
        # bin.loop = 250
        message('Beginning bin loop ',bin.loop)
        
        # Limit plots to -1000ms to +500ms for prod-locked data and -250 to 1000 for stim-locked
        time.boundaries <- 
          list('locked_to_production_onset' = c('start' = -1000, 'end' = 500),
               'locked_to_stimulus_onset' = c('start' = -250, 'end' = 1000))[[lock.loop]]
        times <- times[which((times >= time.boundaries['start']) & (times <= time.boundaries['end']))]
        
        # Define time window bins
        current.windows <- 
          data.frame('start.time' = sort(unique((c(-seq(0, -min(times), by = bin.loop),
                                                   seq(0,max(times) - bin.loop, by = bin.loop))))))
        current.windows$end.time <- current.windows$start.time + bin.loop
        current.windows$start.sample <- time.convert(current.windows$start.time, "times", "samples")
        current.windows$end.sample <- 
          c(time.convert(current.windows$start.time[2:nrow(current.windows)], "times", "samples") - 1,
            current.windows$start.sample[nrow(current.windows)] + time.convert(bin.loop, "times", "samples") - 1)
        current.windows$start.time.label <- 
          gsub("posneg","neg",
               paste0("pos",
                      gsub("-","neg",current.windows$start.time),
                      "ms"))
        current.windows$end.time.label <- 
          gsub("posneg","neg",
               paste0("pos",
                      gsub("-","neg",current.windows$end.time),
                      "ms"))
        current.windows$label <- 
          with(current.windows, 
               paste0(start.time.label, "_to_", end.time.label))
        
        # Initialize data storage
        naming.data.template <- data.frame(matrix(nrow = nrow(trial.info[['pn']]),
                                                  ncol = nrow(current.windows)))
        colnames(naming.data.template) <- current.windows$label
        sentence.data.template <- data.frame(matrix(nrow = nrow(trial.info[['sp']]),
                                                    ncol = nrow(current.windows)))
        colnames(sentence.data.template) <- current.windows$label
        listing.data.template <- data.frame(matrix(nrow = nrow(trial.info[['lp']]),
                                                   ncol = nrow(current.windows)))
        colnames(listing.data.template) <- current.windows$label
        
        
        ## Initialize results storage
        # Template
        data.storage.template <- data.frame(matrix(nrow = length(elec.labels),
                                                   ncol = nrow(current.windows),
                                                   dimnames = list(elec.labels, current.windows$label)))
        data.storage.template$elec <- rownames(data.storage.template)
        data.storage.template$region_clinical <- elec.info[data.storage.template$elec, 'region_clinical']
        data.storage.template$MNI_x <- elec.info[data.storage.template$elec, 'MNI_x']
        data.storage.template$MNI_y <- elec.info[data.storage.template$elec, 'MNI_y']
        data.storage.template$MNI_z <- elec.info[data.storage.template$elec, 'MNI_z']
        
          
        # Differences between tasks
        diff.sentence.naming.ds <- data.storage.template
        diff.sentence.naming.p <- data.storage.template
        diff.sentence.naming.sig <- list()
        diff.sentence.listing.ds <- data.storage.template
        diff.sentence.listing.p <- data.storage.template
        diff.sentence.listing.sig <- list()
        diff.listing.naming.ds <- data.storage.template
        diff.listing.naming.p <- data.storage.template
        diff.listing.naming.sig <- list()
        for(alpha.loop in alphas.to.try){
          diff.sentence.naming.sig[[paste0('alpha=',alpha.loop)]] <- data.storage.template
          diff.sentence.listing.sig[[paste0('alpha=',alpha.loop)]] <- data.storage.template
          diff.listing.naming.sig[[paste0('alpha=',alpha.loop)]] <- data.storage.template
        }; rm(alpha.loop)
        rm(data.storage.template)
        
        ## Loop thru elecs 
        for(elec.loop in 1:length(elec.labels)){ 
          # elec.loop <- 1
          current.elec = elec.labels[elec.loop]
          
          # Set up current.data
          naming.temp.data <- elec.data[['pn']][[current.elec]]
          naming.current.data <- naming.data.template
          sentence.temp.data <- elec.data[['sp']][[current.elec]]
          sentence.current.data <- sentence.data.template
          listing.temp.data <- elec.data[['lp']][[current.elec]]
          listing.current.data <- listing.data.template
          
          for(window.loop in 1:nrow(current.windows)){ 
            # window.loop <- 1
            current.window.start.sample <- current.windows[window.loop,'start.sample']
            current.window.end.sample <- current.windows[window.loop,'end.sample']
            current.window.samples <- current.window.start.sample:current.window.end.sample
            current.window.sample.labels <- time.convert(current.window.samples, "samples", "sample.labels")
            rm(current.window.start.sample, current.window.end.sample, current.window.samples)
            current.window.label <- current.windows[window.loop,'label']
            
            # Fill in data
            naming.current.data[,current.window.label] <- 
              apply(naming.temp.data[,current.window.sample.labels], 1, mean)
            sentence.current.data[,current.window.label] <- 
              apply(sentence.temp.data[,current.window.sample.labels], 1, mean)
            listing.current.data[,current.window.label] <- 
              apply(listing.temp.data[,current.window.sample.labels], 1, mean)
            rm(current.window.sample.labels, current.window.label)
          }; rm(window.loop, naming.temp.data, sentence.temp.data, listing.temp.data)
          
          
          for(sample.loop in colnames(naming.current.data)){ 
            # sample.loop = colnames(naming.current.data)[1]
            
            ## Two-sample Wilcox tests comparing ECoG between tasks
            # Sentence and Naming
            current.wilcox.test <- 
              wilcox.test(sentence.current.data[,sample.loop],
                          naming.current.data[,sample.loop])#, alternative = "greater")
            diff.sentence.naming.p[current.elec, sample.loop] <- current.wilcox.test$p.value
            rm(current.wilcox.test)
            diff.sentence.naming.ds[current.elec, sample.loop] <- 
              cohen.d(sentence.current.data[,sample.loop],
                      naming.current.data[,sample.loop])$estimate
            if(nrow(listing.current.data) > 0){
              # Sentence and Listing
              current.wilcox.test <- 
                wilcox.test(sentence.current.data[,sample.loop],
                            listing.current.data[,sample.loop])#, alternative = "greater")
              diff.sentence.listing.p[current.elec, sample.loop] <- current.wilcox.test$p.value
              rm(current.wilcox.test)
              diff.sentence.listing.ds[current.elec, sample.loop] <- 
                cohen.d(sentence.current.data[,sample.loop],
                        listing.current.data[,sample.loop])$estimate
              # Listing and Naming
              current.wilcox.test <- 
                wilcox.test(listing.current.data[,sample.loop],
                            naming.current.data[,sample.loop])#, alternative = "greater")
              diff.listing.naming.p[current.elec, sample.loop] <- current.wilcox.test$p.value
              rm(current.wilcox.test)
              diff.listing.naming.ds[current.elec, sample.loop] <- 
                cohen.d(listing.current.data[,sample.loop],
                        naming.current.data[,sample.loop])$estimate
            } # if there's any listing data for this patient
            
          }#; rm(sample.loop)
          
          ## Clean up
          rm(current.elec, naming.current.data, sentence.current.data, listing.current.data)
        }; rm(elec.loop)
        
        
        
        ### Correct for multiple comparisons across electrodes
        for(window.loop in current.windows$label){
          # window.loop = current.windows$label[1]
          
          ### P-value correction
          # Sentence vs. listing
          diff.sentence.listing.p[,window.loop] <-
            p.adjust(diff.sentence.listing.p[,window.loop],
                     method = 'fdr')
          
          # Sentence vs. naming
          diff.sentence.naming.p[,window.loop] <-
            p.adjust(diff.sentence.naming.p[,window.loop],
                     method = 'fdr')
          
          # Listing vs. naming
          diff.listing.naming.p[,window.loop] <-
            p.adjust(diff.listing.naming.p[,window.loop],
                     method = 'fdr')
          
        
        ### Get significance
          for(alpha.loop in alphas.to.try){
            # Sentence vs. listing
            diff.sentence.listing.sig[[paste0('alpha=',alpha.loop)]][,window.loop] <-
              as.numeric(diff.sentence.listing.p[,window.loop] < alpha.loop)
            
            # Sentence vs. naming
            diff.sentence.naming.sig[[paste0('alpha=',alpha.loop)]][,window.loop] <-
              as.numeric(diff.sentence.naming.p[,window.loop] < alpha.loop)
            
            # Listing vs. naming
            diff.listing.naming.sig[[paste0('alpha=',alpha.loop)]][,window.loop] <-
              as.numeric(diff.listing.naming.p[,window.loop] < alpha.loop)
          }; rm(alpha.loop)
      }; rm(window.loop)
        
      
        
        
        ###
        ### Save
        ###
        
        current.save.dir <- paste0(path,
                                   'analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/data/',
                                   warp.loop,' data/',
                                   band.loop,
                                   '/patients - individually/',
                                   lock.loop,'/',
                                   bin.loop, 'ms bins/')
        
        ### Cohens' Ds
        # Sentence and listing 
        current.save.diff.dir <- paste0(current.save.dir,'sentence_vs_listing_cohens_d/')
        dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(diff.sentence.listing.ds,
                  paste0(current.save.diff.dir,
                         patient,', ',
                         lock.loop,', sentence_vs_listing_cohens_d.csv'),
                  row.names = FALSE,
                  quote = FALSE)
        # Sentence and naming
        current.save.diff.dir <- paste0(current.save.dir,'sentence_vs_naming_cohens_d/')
        dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(diff.sentence.naming.ds,
                  paste0(current.save.diff.dir,
                         patient,', ',
                         lock.loop,', sentence_vs_naming_cohens_d.csv'),
                  row.names = FALSE,
                  quote = FALSE)
        # Listing and naming
        current.save.diff.dir <- paste0(current.save.dir,'listing_vs_naming_cohens_d/')
        dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(diff.listing.naming.ds,
                  paste0(current.save.diff.dir,
                         patient,', ',
                         lock.loop,', listing_vs_naming_cohens_d.csv'),
                  row.names = FALSE,
                  quote = FALSE)
        
        
        ### ECoG differences between tasks
        for(alpha.loop in alphas.to.try){
          alpha.label = paste0('alpha=',alpha.loop)
          
          # Sentence and listing significant
          current.save.diff.dir <- paste0(current.save.dir,'wilcox_test_sentence_vs_listing_sig_',alpha.label,'/')
          dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
          write.csv(diff.sentence.listing.sig[[alpha.label]],
                    paste0(current.save.diff.dir,
                           patient,', ',
                           lock.loop,', wilcox_test_sentence_vs_listing_sig, ',
                           alpha.label,
                           '.csv'),
                    row.names = FALSE,
                    quote = FALSE)
          
          # Sentence and naming
          current.save.diff.dir <- paste0(current.save.dir,'wilcox_test_sentence_vs_naming_sig_',alpha.label,'/')
          dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
          write.csv(diff.sentence.naming.sig[[alpha.label]],
                    paste0(current.save.diff.dir,
                           patient,', ',
                           lock.loop,', wilcox_test_sentence_vs_naming_sig, ',
                           alpha.label,
                           '.csv'),
                    row.names = FALSE,
                    quote = FALSE)
          
          # Listing and naming
          current.save.diff.dir <- paste0(current.save.dir,'wilcox_test_listing_vs_naming_sig_',alpha.label,'/')
          dir.create(current.save.diff.dir, showWarnings = FALSE, recursive = TRUE)
          write.csv(diff.listing.naming.sig[[alpha.label]],
                    paste0(current.save.diff.dir,
                           patient,', ',
                           lock.loop,', wilcox_test_listing_vs_naming_sig, ',
                           alpha.label,
                           '.csv'),
                    row.names = FALSE,
                    quote = FALSE)
          
        }; rm(alpha.loop)
        
        # Clean up
        rm(save.ecog.dir, save.ecog.t.dir, save.ecog.sig.dir, save.encoding.imp.dir, save.encoding.sig.dir, current.save.which.word.dir, current.save.bayes.dir)
      }#; rm(bin.loop)
      
      # Don't return anything from foreach loop so as not to bog down memory
      return(NULL)
      
    } # patient, lock.loop, warp.loop
  # End parallel processing
  stopCluster(cl)

  beep()  
} # band.loop

message('Script complete. ', Sys.time())

