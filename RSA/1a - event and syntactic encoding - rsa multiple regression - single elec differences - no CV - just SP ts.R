### Do multiple regression on RSA values
### June 2023
### adam.milton.morgan@gmail.com

###
### Readme
###


###
### Setup
###

### Packages
library('beepr') # for playing beeps
library('dplyr') # for data organization, including bind_rows()
library('parallel')
library('foreach')
library('pals')
library('caret')
library('doParallel')
# library('lme4')

# Clean up
rm(list=ls())
cat("\014")
message("Begin multiple regressions on electrode RSA values (correlations across short time windows) on warped data. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 4
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 6
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 14
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 10 # big RDMs
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
# Subset to just good picture naming trials
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add colored text to plots
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Create sigmoids
source(paste0(path,'/analysis/R/functions/sigmoid.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))


### Set up elec info
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
               (elec.info$bad_elec == 0))# &
               # (elec.info$visual_elec == 0) &
               # (elec.info$active == 1))
rownames(elec.info) <- elec.info$patient_elec
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec

keep.stuff <- c(ls(), 'keep.stuff', 'band.loop')

### Loop thru beta/high gamma data
for(band.loop in c('high_gamma','beta')){ ## UNCOMMENT
  # band.loop = c('high_gamma','beta')[2] ## UNCOMMENT
  
  # Clean up
  rm(list = ls()[which(! ls() %in% keep.stuff)])
  
  # Save directory
  output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/')
  
  
  ### Get electrode importances!
  
  ## If already run, skip
  elec.ts.path <- paste0(output.path,'data/1a - elec ts/')
  elec.ts.file <- 'elec RSA multiple regression t-values for PN and SP.RData'
  if(file.exists(paste0(elec.ts.path, elec.ts.file))){
    load(paste0(elec.ts.path, elec.ts.file))
  }else{
    
    # Load data
    message("Loading data... ",Sys.time())
    load(paste0(path,
                'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',band.loop,'/data/elec data with RTs stretched to global median by task - bad trials and .025 to .95 RT outliers excluded.RData')) # loads "stretch.data", "median.rt.samples", "patient.trial.info", and "stretch.sample.labels"
    message('...done! ', Sys.time())
    
    ### Clean up and get metadata
    sp.data <- lapply(stretch.data, function(x){x[['sp']]})
    rm(stretch.data)
    sp.info <- lapply(patient.trial.info, function(x){x[['sp']]})
    rm(patient.trial.info)
    
    ## Get rid of extremely early/late time samples (< 400ms before stim onset and > 500ms post articulation)
    # Get "keep" sample ranges
    sp.keep.sample.range <- c((-median.rt.samples['sp']) - time.convert(400, "times", "samples"),
                              time.convert(500, "times", "samples"))
    # Subset
    sp.data <- lapply(sp.data, function(x){
      lapply(x, function(y){
        current.samples <- time.convert(colnames(y), "sample.labels", "samples")
        keep.cols <- which((current.samples >= sp.keep.sample.range[1]) & (current.samples <= sp.keep.sample.range[2]))
        return(y[,keep.cols])
      })
    })
    
    # Get new sample labels
    sp.sample.labels <- colnames(sp.data[[1]][[1]])
    
    # Subset data to just patients who produced passives
    patients.who.produced.passives <- sapply(sp.info, function(x){sum(x$verb_voice == 'passive') > 4})
    patients.who.produced.passives <- names(patients.who.produced.passives[patients.who.produced.passives])
    sp.data <- sp.data[patients.who.produced.passives]
    
    # Good elecs:
    patient.elecs <- lapply(sp.data, names)
    patient.elecs <- lapply(patient.elecs, function(x){x[x %in% use.these.elecs]})
    elecs.to.loop.thru <- data.frame('elec' = unlist(patient.elecs))
    elecs.to.loop.thru$patient <- substr(elecs.to.loop.thru$elec, 1, 5)
    rownames(elecs.to.loop.thru) <- elecs.to.loop.thru$elec
    elecs.to.loop.thru$region_clinical <- elec.info[elecs.to.loop.thru$elec, 'region_clinical']
    
    # Subset data to just use.these.elecs
    for(patient in names(sp.data)){
      sp.data[[patient]] <- sp.data[[patient]][patient.elecs[[patient]]]
    }; rm(patient)
    gc()
    
    # Load stimulus image info
    stimulus.info <- read.csv(paste0(path, 'stimulus creation/stimulus_info_lookup_for_analysis.csv'))[,c('item','verb','active.subject','active.object','image.id.active.subject.left','image.id.active.subject.right','finite.active.verb','finite.passive.verb')]
    
    # Flesh out scene info
    sp.info <- lapply(sp.info, function(x){
      x$active.subject.position <- ifelse(grepl("right", x$image_id), "right", "left")
      x$event.semantics <- gsub("-left", "", gsub("-right", "", x$image_id, fixed = TRUE), fixed = TRUE)
      return(x)
    })
    
    
    # ### Add event semantic features
    # # Load event semantic features
    # # event.categories <- read.csv(paste0(path, 'analysis/R/define event semantic features/GLoVe word embeddings/output/event_semantic_features.csv'))
    # event.categories <- read.csv(paste0(path, 'analysis/R/define event semantic features/GLoVe word embeddings/output/event_semantic_RDM_GloVe_wikipedia.csv'))
    # rownames(event.categories) <- event.categories$verb1_verb2
    # 
    # ## Subset to just trials with good verbs
    # # SP
    # sp.patient.keep.trials <- lapply(sp.info, function(x){which(x$verb %in% event.categories$verb1)})
    # for(patient in names(sp.info)){
    #   if(! is.null(sp.data[[patient]])){
    #     sp.data[[patient]] <- lapply(sp.data[[patient]], function(x){
    #       droplevels(x[sp.patient.keep.trials[[patient]],])
    #     })  
    #   }
    #   sp.info[[patient]] <- droplevels(sp.info[[patient]][sp.patient.keep.trials[[patient]],])
    # }; rm(patient, sp.patient.keep.trials)
    
    
    ### Load GPT-2 sentence (contextual) embeddings (vectors)
    contextual.embeddings <- read.csv(paste0(path, 'analysis/R/define event semantic features/GPT2 contextual embeddings/output/stimuli_contextual_embeddings_from_GPT2_layer8.csv'))
    colnames(contextual.embeddings) <- gsub(".", "-", colnames(contextual.embeddings), fixed = TRUE)
    # Get correlations
    event.cors <- cor(contextual.embeddings)
    # scale -1 to 1
    event.cors <- event.cors - min(event.cors)
    event.cors <- ((event.cors / max(event.cors)) * 2) - 1
    event.cors <- 
      data.frame('event1' = rownames(event.cors)[row(event.cors)], 
                 'event2' = colnames(event.cors)[col(event.cors)], 
                 'cor' = c(event.cors))
    event.cors$diff <- 1 - event.cors$cor
    event.cors$diff <- event.cors$diff - 1 # center... so really just negative (stretched) correlation
    rownames(event.cors) <- paste0(event.cors$event1, "_", event.cors$event2)
    
    
    ## Subset to just trials with good verbs
    sp.patient.keep.trials <- lapply(sp.info, function(x){which(x$event.semantics %in% event.cors$event1)})
    for(patient in names(sp.info)){
      if(! is.null(sp.data[[patient]])){
        sp.data[[patient]] <- lapply(sp.data[[patient]], function(x){
          droplevels(x[sp.patient.keep.trials[[patient]],])
        })
      }
      sp.info[[patient]] <- droplevels(sp.info[[patient]][sp.patient.keep.trials[[patient]],])
    }; rm(patient, sp.patient.keep.trials)
    
    
    ### Smoothing hyperparameters
    half.n.smoothing.window.ms = 100
    samples <- time.convert(sp.sample.labels, "sample.labels", "samples")
    sample.diffs <- c()
    for(i in 2:length(samples)){sample.diffs <- c(sample.diffs, samples[i] - samples[i-1])}
    half.n.smoothing.samples <- time.convert(half.n.smoothing.window.ms, "times", "samples") / mean(sample.diffs)
    
    
    ###
    ### Get t-values
    ###
    
    ### Loop thru patients, elecs, trials
    message('Beginning cor calculations: ',
            ' running groups in parallel with ',n.cores.to.use,' cores (',band.loop,'). ',
            Sys.time())
    
    
    ## Set up parallel processing
    # Close any old parallel backends
    unregister_dopar()
    # Set up parallel workers in case caret wants to use it
    cl <- makeCluster(n.cores.to.use, type = "FORK")
    registerDoParallel(cl)
    
    term.stats <- 
      foreach(elec.loop = elecs.to.loop.thru$elec) %dopar% {
        # elec.loop = elecs.to.loop.thru$elec[1]
        
        # Set seed
        current.seed = 10000 + which(elecs.to.loop.thru$elec == elec.loop)
        set.seed(seed = current.seed)
        
        ## Set up
        patient <- elecs.to.loop.thru[elec.loop, 'patient']
        
        # Get the indices of a hypothetical difference matrix
        indices.sp <- which(lower.tri(matrix(nrow = nrow(sp.info[[patient]]),
                                             ncol = nrow(sp.info[[patient]]))), arr.ind = TRUE)
        
        ### Get elec distances for Sentence Production
        # Loop thru time windows
        sp.elec.distances <- list()
        for(sample.loop in sp.sample.labels){
          # sample.loop = sp.sample.labels[1]
          
          # Initialize if first data for this elec
          if(length(sp.elec.distances) == 0){
            sp.elec.distances[['metadata']] <- 
              data.frame('voice1' = sp.info[[patient]]$verb_voice[indices.sp[,'row']],
                         'voice2' = sp.info[[patient]]$verb_voice[indices.sp[,'col']],
                         'event1' = sp.info[[patient]]$event.semantics[indices.sp[,'row']],
                         'event2' = sp.info[[patient]]$event.semantics[indices.sp[,'col']],
                         'word1' = sp.info[[patient]]$word[indices.sp[,'row']],
                         'word2' = sp.info[[patient]]$word[indices.sp[,'col']],
                         'rt1' = sp.info[[patient]]$production_latency[indices.sp[,'row']],
                         'rt2' = sp.info[[patient]]$production_latency[indices.sp[,'col']])
          } # if(length(sp.elec.distances) == 0){
          
          # Turn distances into dataframes
          sp.elec.distances[[sample.loop]] <- 
            data.frame('temp' = abs(sp.data[[patient]][[elec.loop]][,sample.loop][indices.sp[,'row']] - 
                                      sp.data[[patient]][[elec.loop]][,sample.loop][indices.sp[,'col']]))
          
          # Rename column and add to data
          colnames(sp.elec.distances[[sample.loop]]) <- sample.loop
          
        }; rm(sample.loop)
        
        rm(indices.sp, patient)
        
        # Combine into one dataframe
        sp.elec.distances <- bind_cols(sp.elec.distances)
        
        
        ### Create IVs
        sp.ivs <- data.frame(
          'diff.voice' = as.numeric(sp.elec.distances$voice1 != sp.elec.distances$voice2),
          'diff.1st.word'= as.numeric(sp.elec.distances$word1 != sp.elec.distances$word2),
          'diff.event.semantics' = event.cors[paste0(sp.elec.distances$event1,"_",sp.elec.distances$event2), "diff"],
          'diff.log.rt' = abs(log(sp.elec.distances$rt1) - log(sp.elec.distances$rt2)))
        # Define "model.terms"
        sp.model.terms <- colnames(sp.ivs)
        
        # Add IVs to data
        sp.elec.distances <- cbind(sp.ivs, sp.elec.distances)
        
        
        ### Stats
        ## Sentence Production
        # Set up storage
        sp.ts <- data.frame(matrix(nrow = length(sp.sample.labels),
                                   ncol = length(sp.model.terms)))
        colnames(sp.ts) <- sp.model.terms
        rownames(sp.ts) <- sp.sample.labels
        
        # Get stats for each window
        for(sample.loop in sp.sample.labels){
          # sample.loop = sp.sample.labels[1]
          
          # Set up data
          current.data <- sp.elec.distances[, c(sample.loop, sp.model.terms)]  
          names(current.data)[which(names(current.data) == sample.loop)] <- 'distance'
          
          # Evaluate model.terms one at a time
          model <- summary(lm(distance ~ ., data = current.data))$coefficients
          
          # Store model results
          for(term.loop in rownames(model)[rownames(model) != "(Intercept)"]){
            sp.ts[sample.loop, term.loop] <- model[term.loop, 't value']
          }; rm(term.loop)
          rm(current.data, model)
          
        }; rm(sample.loop)
        rm(sp.elec.distances)
        
        
        ### Smooth
        sp.ts.smooth <- apply(sp.ts, 2, function(x){return(smoothing(x, n.samples.pre = half.n.smoothing.samples))})
        # Remove NAs introduced by smoothing
        sp.ts.smooth <- sp.ts.smooth[-c(1:half.n.smoothing.samples,
                                        (nrow(sp.ts.smooth) - half.n.smoothing.samples + 1):nrow(sp.ts.smooth)),]
        
        ## Output
        return.stats <- list('sp.ts' = sp.ts,
                             'sp.ts.smooth' = sp.ts.smooth,
                             'elec' = elec.loop)
        return(return.stats)
        
      } # elec.loop
    
    ## Close parallel backend
    stopCluster(cl)
    unregister_dopar()
    
    message("Done w stats! ",Sys.time())
    beep()
    
    
    ## Reorganize
    # Label output with elec names
    for(elec.loop in 1:length(term.stats)){
      # elec.loop = 1
      names(term.stats)[elec.loop] <- term.stats[[elec.loop]]$elec
    }; rm(elec.loop)
    
    # Separate t-values and p-values
    sp.ts <- lapply(term.stats, function(x){return(x[['sp.ts']])})
    sp.ts.smooth <- lapply(term.stats, function(x){return(x[['sp.ts.smooth']])})
    
    # Save term names
    sp.model.terms <- colnames(sp.ts[[1]])
    
    
    ### Save
    save.these <- c(
      'sp.ts',
      'sp.ts.smooth',
      'sp.sample.labels',
      'sp.model.terms',
      'median.rt.samples'
    )
    
    dir.create(elec.ts.path, showWarnings = FALSE, recursive = TRUE)
    save(list = save.these,
         file = paste0(elec.ts.path, elec.ts.file))
    rm(save.these, sp.data, sp.info, term.stats)
    
    
    ### Save smoothing hps
    save.these.smoothing <- c(
      'half.n.smoothing.samples',
      'half.n.smoothing.window.ms'
    )
    
    save(list = save.these.smoothing,
         file = paste0(elec.ts.path, 'RSA smoothing hyperparameters.RData'))
    rm(save.these.smoothing)
    
    
  } # If t-values already run, skip
  
}; rm(band.loop)