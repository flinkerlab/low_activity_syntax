### Cluster RSIs (zs.adjusted) trained on shuffled data to determine how much signal NMF can find in noise
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
library('NMF') # for NMF
# library('lme4')

# Clean up
rm(list=ls())
# cat("\014")
message("Begin multiple regressions on electrode RSA values (correlations across short time windows) on warped data. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 4
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 7
    n.cores.to.use.nmf = 7
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 15
    n.cores.to.use.nmf = 15
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 10 # big RDMs
  n.cores.to.use.nmf = 32
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

## Define function for squishing data: reduce extreme low and high vals
squish <- function(
    vals.to.squish, 
    scaling.factor = qnorm(.01, lower.tail = FALSE)
){
  # Get indices of NAs and temporarily replace them with 1s, to be switched back at end
  na.indices <- which(is.na(vals.to.squish))
  vals.to.squish[na.indices] <- 1
  # Scale and squish and unscale
  out <- vals.to.squish / scaling.factor # Scale so significant values > 1, n.s. < 1
  out[out > 1] <- log(out[out > 1]) + 1 # Reduce extremely high values
  out[out < 1] <- (exp(out[out < 1]) / exp(1)) ^ 3 # Reduce non-significant values
  out <- out * scaling.factor # Unscale  
  out[na.indices] <- NA
  return(out)
}


### Loops!
## Clean up
keep.all.this1 <- c(ls(), 
                   'keep.all.this1', 
                   'band.loop')

## Former loops
band.loop = c('high_gamma','beta')[1]

 rm(list = ls()[! ls() %in% keep.all.this1])
    gc()
    
   
## Directories
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/', band.loop, '/')
elec.shuffle.zs.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/shuffle zs and sig elecs for each shuffle loop/')
elec.shuffle.zs.files <- list.files(elec.shuffle.zs.path)

## Shuffle loops
n.shuffle.loops <- length(elec.shuffle.zs.files)

### Load elec means
load(paste0(path,'analysis/R/electrode means/means - warped data - multiband/output - high_gamma/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses
rm(elec.sds, elec.ses)

# Get metadata
tasks <- names(elec.means)

# Add NA dataframe for Patient006 LP block
blank.Patient006.lp.df <- data.frame(matrix(nrow = nrow(elec.means$lp[[1]]),
                                       ncol = ncol(elec.means$pn[['Patient006']])))
colnames(blank.Patient006.lp.df) <- colnames(elec.means$pn[['Patient006']])
rownames(blank.Patient006.lp.df) <- rownames(elec.means$lp[[1]])
elec.means$lp[['Patient006']] <- blank.Patient006.lp.df
rm(blank.Patient006.lp.df)

# Collapse across patients
elec.means <- lapply(elec.means, bind_cols)


### Get HGA part of stacked data
## Only patients that produced passives
attach(paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/adjusted zs with sig windows and elecs.RData'))
elecs.with.passive.data <- names(which(apply(data.frame(bind_cols(lapply(zs.adjusted, function(x){return(x$diff.voice)}))), 2, function(y){!any(is.na(y))})))
detach()
patients.who.produced.passives <- unique(substr(elecs.with.passive.data, 1, 5))

### HGA
  ### Set up elec means data
  # Get sample labels to keep by task
  og.samples <- list()
  keep.these.samples <- list()
  keep.these.sample.labels <- list()
  for(task.loop in tasks){ # task.loop = tasks[1]
    # Get the original column samples
    og.samples[[task.loop]] <- 
      time.convert(rownames(elec.means[[task.loop]]), "sample.labels", "samples")
    # Subset to just those that come after stim onset
    keep.these.samples[[task.loop]] <-
      og.samples[[task.loop]][which(og.samples[[task.loop]] > (-median.rt.samples[task.loop]))]
    # Convert back to sample labels (column labels)
    keep.these.sample.labels[[task.loop]] <- 
      time.convert(keep.these.samples[[task.loop]], "samples", "sample.labels")
  }; rm(task.loop, og.samples, keep.these.samples)
  
  # Smoothing hyperparameters
  load(paste0(output.path, 'data/1a - elec ts/RSA smoothing hyperparameters.RData'))
  
  # Remove pre-stim data
  post.stim.elec.means <- list()
  for(task.loop in tasks){
    # task.loop = tasks[1]
    post.stim.elec.means[[task.loop]] <- 
      elec.means[[task.loop]][keep.these.sample.labels[[task.loop]], elecs.with.passive.data]
    
    # Smooth
    half.n.hga.smoothing.samples <- ceiling(half.n.smoothing.samples / 2)
    for(col.loop in 1:ncol(post.stim.elec.means[[task.loop]])){
      # col.loop = 1
      post.stim.elec.means[[task.loop]][,col.loop] <- 
        smoothing(post.stim.elec.means[[task.loop]][,col.loop],
                  n.samples.pre = half.n.hga.smoothing.samples)
    }; rm(col.loop)
    post.stim.elec.means[[task.loop]] <- 
      post.stim.elec.means[[task.loop]][-c(1:half.n.hga.smoothing.samples,
                                           (nrow(post.stim.elec.means[[task.loop]]) - half.n.hga.smoothing.samples + 1):nrow(post.stim.elec.means[[task.loop]])),]
  }; rm(task.loop)
  
  ### Stack
  hga.stacked.data <- bind_rows(post.stim.elec.means)
  rm(post.stim.elec.means)
  
  ## Do the same transformation you did to the t-values
  hga.stacked.data <- data.frame(apply(hga.stacked.data, 2, function(x){squish(x, 1)}))


### Clean up
keep.all.this2 <- c(ls(), 
                   'keep.all.this2', 
                   'alpha.hps',
                   'alpha.loop')


### Set up directories
# Significance hyperparameters to loop through
alpha.hps <- c('including_inactive_only_if_RSA_sig_alpha=.05_min=100ms',
               'including_inactive_only_if_RSA_sig_alpha=.01_min=50ms')

# for(alpha.loop in alpha.hps[2]){
alpha.loop = alpha.hps[1]

rm(list = ls()[! ls() %in% keep.all.this2])
gc()

# Alpha threshold
  alpha.threshold <- gsub("including_inactive_only_if_RSA_sig_","",alpha.loop)


# Save dir
shuffle.nmf.path <- paste0(output.path,'data/1f - shuffled data nmf results/',
                           'with depth elecs','/',
                           'with inactive elecs','/',
                           "with hga/",
                           alpha.loop,'/',
                           'model weights/')

# How many shuffle loops already done?
shuffle.loops.completed <- 
  as.numeric(gsub("nmf - shuffled data - loop ","",
                  gsub(".RData","",list.files(shuffle.nmf.path), fixed = TRUE), fixed = TRUE))
shuffle.loops.to.do <- (1:n.shuffle.loops)[which(! (1:n.shuffle.loops) %in% shuffle.loops.completed)]

### Set up elec info
elec.info <- 
  read.csv(paste0(path,
                  '/analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
elec.info$depth_elec <- as.numeric(substr(elec.info$elec,1,1) == "D")
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


# Define which elecs to use
  use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec  


###
### Stack data and NMF
###

### Loop thru patients, elecs, trials
message('Beginning cor calculations: ',
        ' running groups in parallel with ',n.cores.to.use,' cores (',band.loop,'). ',
        Sys.time())

### Run shuffle loops!
for(shuffle.loop in shuffle.loops.to.do){
  # shuffle.loop = shuffle.loops.to.do[1]
  start.time <- Sys.time()
  
  set.seed(seed = shuffle.loop * 1000 + 401)
  
  # Load shuffled data
  load(paste0(elec.shuffle.zs.path, elec.shuffle.zs.files[shuffle.loop])) # loads "shuffle.zs.adjusted" and "shuffle.zs.sig.elecs"
  
  ##
  ## Stack data
  ##
  
  ### Add terms ts from SP
  ## Get sample labels to keep by task
  sp.og.samples <- time.convert(rownames(shuffle.zs.adjusted[[1]]), "sample.labels", "samples")
  sp.keep.these.samples <- sp.og.samples[which(sp.og.samples > (-median.rt.samples['sp']))]
  sp.keep.these.sample.labels <- time.convert(sp.keep.these.samples, "samples", "sample.labels")
  
  ## Get SP data to stack
  # diff.voice
  sp.voice <- data.frame(bind_cols(lapply(shuffle.zs.adjusted, function(x){
    x[sp.keep.these.sample.labels, 'diff.voice']
  })))
  sp.voice <- sp.voice[,colnames(sp.voice)[colnames(sp.voice) %in% elecs.with.passive.data]]
  # diff.1st.word
  sp.1st.word <- data.frame(bind_cols(lapply(shuffle.zs.adjusted, function(x){
    x[sp.keep.these.sample.labels, 'diff.1st.word']
  })))
  sp.1st.word <- sp.1st.word[,colnames(sp.1st.word)[colnames(sp.1st.word) %in% elecs.with.passive.data]]
  # diff.event.semantics
  sp.event.semantics <- data.frame(bind_cols(lapply(shuffle.zs.adjusted, function(x){
    x[sp.keep.these.sample.labels, 'diff.event.semantics']
  })))
  sp.event.semantics <- sp.event.semantics[,colnames(sp.event.semantics)[colnames(sp.event.semantics) %in% elecs.with.passive.data]]
  
  ## Stack
  stacked.data <- bind_rows(list(sp.voice,
                                 sp.1st.word,
                                 sp.event.semantics))
  rm(sp.og.samples, sp.keep.these.samples, sp.keep.these.sample.labels, sp.voice, sp.1st.word, sp.event.semantics)
  
  
  ### Combine HGA and z's
    # Subset to just electrodes you're using
    elecs.to.use <- names(shuffle.zs.adjusted)[names(shuffle.zs.adjusted) %in% elecs.with.passive.data]
    hga.stacked.data <- hga.stacked.data[,elecs.to.use]
    
    # Combine
    stacked.data <- bind_rows(list(stacked.data, hga.stacked.data))
  
  
  ### Impute missing data (Patient006's LP block and passive term for Patient004)
  elecs.w.missing.data <- names(which(apply(stacked.data, 2, function(x){any(is.na(x))})))
  for(elec.with.nas.loop in elecs.w.missing.data){
    # elec.with.nas.loop = elecs.w.missing.data[1]
    missing.rows <- which(is.na(stacked.data[,elec.with.nas.loop]))
    elecs.with.those.rows <- names(which(apply(stacked.data[missing.rows,], 2, function(x){! any(is.na(x))})))
    elec.rs <- apply(stacked.data[-missing.rows, elecs.with.those.rows], 2, function(z){cor(x = z, y = stacked.data[-missing.rows, elec.with.nas.loop])})
    most.similar.elec <- names(elec.rs[which.max(elec.rs)])
    
    # Replace NAs
    stacked.data[missing.rows, elec.with.nas.loop] <-
      stacked.data[missing.rows, most.similar.elec]
    
    # Clean up
    rm(elec.rs, missing.rows, elecs.with.those.rows, most.similar.elec)
  }; rm(elec.with.nas.loop, elecs.w.missing.data)
  
  
  ## Select appropriate elecs depending on loops
  if(alpha.loop == "including_inactive_only_if_RSA_sig_alpha=.05_min=100ms"){
    all.active.elecs <- elec.info[(elec.info$use.these.elecs == 1) &
                                    (elec.info$active == 1),]$patient_elec
    all.sig.elecs <- unique(unlist(shuffle.zs.sig.elecs[['alpha=.05_min=100ms']], use.names = FALSE))
    sp.sig.elecs.all <- intersect(names(shuffle.zs.adjusted), union(all.active.elecs, all.sig.elecs))
  }
  if(alpha.loop == "including_inactive_only_if_RSA_sig_alpha=.01_min=50ms"){
    all.active.elecs <- elec.info[(elec.info$use.these.elecs == 1) &
                                    (elec.info$active == 1),]$patient_elec
    all.sig.elecs <- unique(unlist(shuffle.zs.sig.elecs[['alpha=.01_min=50ms']], use.names = FALSE))
    sp.sig.elecs.all <- intersect(names(shuffle.zs.adjusted), union(all.active.elecs, all.sig.elecs))
  }
  stacked.data <- stacked.data[, colnames(stacked.data)[colnames(stacked.data) %in% sp.sig.elecs.all]]
  
  
  message("Running NMF models! ", Sys.time())
  # Run NMF
  nmf.ranks.to.try <- 5
  nmf.models <- nmf(stacked.data, 
                    rank = nmf.ranks.to.try, 
                    nrun = n.cores.to.use.nmf, # high - default is 30 (i.e., runs NMF 30 times per rank to find a consensus matrix)
                    .opt = 'vP', # mac (already running in parallel automatically)
                    seed = 123)
  
  # Clean up
  rm(stacked.data)
  
  
  
  ### Loop thru NMFS and create data for brain plots
  nmf.rank.loop <- nmf.ranks.to.try
  # nmf.rank.loop = nmf.ranks.to.try
  current.label = paste0('rank=',nmf.rank.loop)  
  
  # Get data from results
  shuffle.weights <- data.frame(t(coef(nmf.models)))#$fit[[as.character(nmf.rank.loop)]])))
  colnames(shuffle.weights) <- paste0("NMF_",1:ncol(shuffle.weights))
  
  # Assign elecs to clusters corresponding to the factors - whichever factor has highest value
  cluster.assignments <- data.frame('electrode' = row.names(shuffle.weights),
                                    'cluster' = apply(shuffle.weights, 1, which.max))
  # cluster.assignments[names(which(apply(shuffle.weights, 1, max) < cluster.threshold)),]$cluster <- nmf.rank.loop + 1
  
  # Combine
  colnames(shuffle.weights) <- paste0('weight_',colnames(shuffle.weights))
  shuffle.weights <- cbind(cluster.assignments, shuffle.weights)
  
  
  ###
  ### Save
  ###
  save.these <- c('shuffle.weights')
  
  dir.create(shuffle.nmf.path, showWarnings = FALSE, recursive = TRUE)
  save(list = save.these,
       file = paste0(shuffle.nmf.path, 'nmf - shuffled data - loop ',shuffle.loop,'.RData'))
  rm(save.these, shuffle.zs.adjusted, shuffle.zs.sig.elecs, shuffle.weights, cluster.assignments, nmf.models)
  gc()
  
  
  print(paste0("Significance loop ",alpha.loop," (",which(alpha.hps == alpha.loop)," of ",length(alpha.hps),"); ","Shuffle loop ",shuffle.loop," (",which(shuffle.loops.to.do == shuffle.loop)," of ",length(shuffle.loops.to.do),") completed. Duration: ", round(Sys.time() - start.time, 2)))
  
  
}; rm(shuffle.loop)
# }; rm(alpha.loop)

# }; rm(band.loop)



# Finish!
message('Script completed successfully. ',Sys.time())







