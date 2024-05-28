### Shuffle elec-cluster assignments and compare actual cluster means to shuffled cluster means, original RSA model
### October 2023
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
library('Hmisc') # for wtd.mean and wtd.var
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
    n.cores.to.use.nmf = 8
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 16
    n.cores.to.use.nmf = 20
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
# Get min and max
source(paste0(path,'/analysis/R/functions/min_max.R'))

### Load plotting colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

### Set up 
## Elec info
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
               (elec.info$bad_elec == 0) &
               # (elec.info$visual_elec == 0) &
               (elec.info$active == 1))
rownames(elec.info) <- elec.info$patient_elec
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec

### Clean up
keep <- c(ls(), 
          'band.loop', 
          'alpha.hps',
          'sig.hp.loop')

### Loop thru stuff
# for(band.loop in c('high_gamma','beta')){ ## UNCOMMENT
band.loop = c('high_gamma','beta')[1] ## UNCOMMENT


# Define alpha hyperparameters to loop thru
  alpha.hps <- c('including_inactive_only_if_RSA_sig_alpha=.05_min=100ms',
                 'including_inactive_only_if_RSA_sig_alpha=.01_min=50ms')

# for(sig.hp.loop in alpha.hps){
sig.hp.loop = 'including_inactive_only_if_RSA_sig_alpha=.05_min=100ms'

# Clean up
rm(list = ls()[which(!ls() %in% keep)])
gc()

# Directories
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/', band.loop, '/')
include.depth.dir <- 'with depth elecs'
include.inactive.dir <- 'with inactive elecs'
  
  
  ###
  ### Define feature spaces
  ###
  
  groupings <- list()
  weights <- list()
  
  ##
  ## Add regions:
  ##
  
  attach(paste0(path,'analysis/R/define ROIs/output/all region definitions.RData'))
  gross.rois <- gross.rois
  detach()
  
  # Limit
  gross.rois <- gross.rois[c('STG','MTG','pericentral','IFG','MFG','IPL','visual')]
  groupings[['rois']] <- list()
  for(roi.loop in names(gross.rois)){
    groupings[['rois']][[roi.loop]] <-
      elec.info[(elec.info$patient_elec %in% use.these.elecs) &
                  (elec.info$region_clinical %in% gross.rois[[roi.loop]]), 'patient_elec']
  }; rm(roi.loop)
  
  
  ##
  ## NMF clusters
  ##
  
  
  ## Loop thru sig elecs
  for(sig.elec.hp.loop in alpha.hps[1:2]){
    # sig.elec.hp.loop = alpha.hps[1]
  
  nmf.results.path <- paste0(output.path, 'data/1e - nmf/',
                             include.depth.dir,'/',
                             include.inactive.dir,'/',
                             "with hga/",
                             sig.elec.hp.loop,'/',
                             '/individual models/csvs for brain plots/')
  nmf.ranks <- sort(as.numeric(gsub("rank=","",list.files(nmf.results.path))))
  
  # Loop thru ranks
  # for(rank.loop in nmf.ranks){
  rank.loop = 5
  rank.label = paste0('rank=',rank.loop,' - sig elec threshold=',sig.elec.hp.loop)
  
  # Get cluster elecs
  groupings[[rank.label]] <- list()
  
  # Load cluster assignments
  cluster.assignments <- read.csv(paste0(nmf.results.path, 'rank=', rank.loop, '/cluster_assignments.csv'))
  
  # Get cluster weight for each electrode 
  weights[[rank.label]] <- read.csv(paste0(nmf.results.path, 'rank=', rank.loop, '/weights.csv'))
  weights[[rank.label]]$weight <- NA
  for(i in 1:nrow(weights[[rank.label]])){
    # i = 1
    weights[[rank.label]]$weight[i] <- weights[[rank.label]][i, which.max(unlist(weights[[rank.label]][i,1:rank.loop]))]
  }; rm(i)
  weights[[rank.label]] <- cbind(cluster.assignments[, c('electrode', 'cluster')], weights[[rank.label]][, 'weight', drop = FALSE])
  row.names(weights[[rank.label]]) <- weights[[rank.label]]$electrode
  
  # Split cluster assignments by cluster
  cluster.assignments <- split(cluster.assignments, cluster.assignments$cluster)
  
  # Get elecs in each cluster
  for(cluster.loop in names(cluster.assignments)){
    groupings[[rank.label]][[cluster.loop]] <- cluster.assignments[[cluster.loop]]$electrode
  }; rm(cluster.loop)
  
  
  # }; rm(rank.loop)
  }; rm(sig.elec.hp.loop)
  
  
  ### Re-order
  # groupings <- groupings[c(4,1,5,3,2,6:length(groupings))]
  
  
  ###
  ### Load t-values
  ###
  
  ### Get electrode importances!
  # # t-values
  # load(paste0(output.path,'data/elec ts/elec RSA multiple regression t-values for PN and SP.RData'))
  # Adjusted t-values
  load(paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/adjusted zs with sig windows and elecs.RData')) # loads significance.hps, zs.adjusted, zs.sig.elecs, zs.sig.windows
  terms <- names(zs.adjusted[[1]])
  rm(significance.hps, zs.sig.elecs, zs.sig.windows)
  
  # RTs
  load(paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
  median.rt.times <- time.convert(median.rt.samples, "samples", "times")
  names(median.rt.times) <- names(median.rt.samples)
  
  
  ### Loop thru groupings
  for(grouping.loop in names(groupings)){ ## UNCOMMENT
    # (grouping.loop = names(groupings)[1]) ## UNCOMMENT
    
    # Update
    message('Calculating significance and plotting ',grouping.loop,' (grouping ',which(names(groupings) == grouping.loop),' of ',length(groupings),'). ',Sys.time())
    
    # All elecs in this grouping loop
    all.elecs.this.grouping.loop <- unique(unlist(groupings[[grouping.loop]], recursive = TRUE, use.names = FALSE))
    all.patients.this.grouping.loop <- unique(substr(all.elecs.this.grouping.loop, 1, 5))
    all.groups.this.grouping.loop <- names(groupings[[grouping.loop]])
    
    ### Get means and SEs
    sp.means <- list()
    sp.ses <- list()
    sp.means.wtd <- list()
    sp.ses.wtd <- list()
    
    ## Loop through groups (clusters/ROIs and get means/SES)
    for(group.loop in names(groupings[[grouping.loop]])){
      # group.loop = names(groupings[[grouping.loop]])[1]
      
      # Current elecs
      current.elecs <- groupings[[grouping.loop]][[group.loop]]
      current.elecs <- current.elecs[current.elecs %in% names(zs.adjusted)]
      if(! is.null(weights[[grouping.loop]])){ # If weights
        current.weights <- weights[[grouping.loop]][current.elecs, 'weight']
      } # If weights
      
      ### Get mean t-values per cluster, SP
      sp.means[[group.loop]] <- list()
      sp.ses[[group.loop]] <- list()
      sp.means.wtd[[group.loop]] <- list()
      sp.ses.wtd[[group.loop]] <- list()
      for(sp.loop in terms){
        # sp.loop = terms[1]
        
        # Get data
        current.data <- do.call(cbind, lapply(zs.adjusted[current.elecs], function(x){return(x[,sp.loop,drop = FALSE])}))
        
        # Get stats
        sp.means[[group.loop]][[sp.loop]] <- apply(current.data, 1, mean, na.rm = TRUE)
        sp.ses[[group.loop]][[sp.loop]] <- apply(current.data, 1, sd) / sqrt(ncol(current.data))
        if(! is.null(weights[[grouping.loop]])){ # If weights
          sp.means.wtd[[group.loop]][[sp.loop]] <- apply(current.data, 1, wtd.mean, weights = current.weights, normwt = TRUE, na.rm = TRUE)  
          sp.ses.wtd[[group.loop]][[sp.loop]] <- apply(current.data, 1, wtd.var, weights = current.weights, normwt = TRUE, na.rm = TRUE) / sqrt(ncol(current.data))
        } # If weights
        
        rm(current.data)
      }; rm(sp.loop)
    }; rm(group.loop)
    
    
    ###
    ### Get shuffle distributions
    ###
    
    # Storage
    sp.shuff.data <- list()
    sp.shuff.data.wtd <- list()
    
    # How many shuffle loops?
    n.shuffles <- 1000
    
    ## Loop through groups (clusters/ROIs and get means/SES)
    for(group.loop in names(groupings[[grouping.loop]])){
      # group.loop = names(groupings[[grouping.loop]])[1]
      
      message("Getting shuffle means for ",group.loop," (group ",which(names(groupings[[grouping.loop]]) == group.loop)," of ",length(groupings[[grouping.loop]]),"). ",Sys.time())
      
      ## Set up parallel processing
      # Close any old parallel backends
      unregister_dopar()
      # Set up parallel workers in case caret wants to use it
      cl <- makeCluster(n.cores.to.use, type = "FORK")
      registerDoParallel(cl)
      
      ## Get shuff means in parallel
      shuff.means <- 
        foreach(shuffle.loop = 1:n.shuffles,
                .options.nws = list(chunkSize = 5)) %dopar% {
                  # shuffle.loop = 1
                  
                  # Seed
                  set.seed(seed = 1000*shuffle.loop)
                  
                  # Shuffle groupings
                  shuff.assignments <- stack(groupings[[grouping.loop]])
                  colnames(shuff.assignments) <- c('elec','group')
                  shuff.assignments$group <- sample(shuff.assignments$group)
                  shuff.assignments <- split(shuff.assignments$elec, shuff.assignments$group)
                  
                  # Current elecs
                  current.elecs <- shuff.assignments[[group.loop]]
                  current.elecs <- current.elecs[current.elecs %in% names(zs.adjusted)]
                  
                  ### Get mean t-values per cluster
                  sp.current.shuff.means <- list()
                  sp.current.shuff.means.wtd <- list()
                  for(sp.term.loop in terms){
                    current.data <- do.call(cbind, lapply(zs.adjusted[current.elecs], function(x){return(x[,sp.term.loop,drop = FALSE])}))
                    sp.current.shuff.means[[sp.term.loop]] <- apply(current.data, 1, mean, na.rm = TRUE)
                    
                    # Weighted mean
                    if(! is.null(weights[[grouping.loop]])){ # If weights
                      current.weights <- weights[[grouping.loop]][current.elecs, 'weight']
                      sp.current.shuff.means.wtd[[sp.term.loop]] <- apply(current.data, 1, wtd.mean, weights = current.weights, normwt = TRUE, na.rm = TRUE)  
                    } # If weights
                    
                    # Clean up
                    rm(current.data)
                  }; rm(sp.term.loop)
                  
                  # Clean up
                  rm(shuff.assignments, current.elecs)
                  
                  # Output
                  return(list('sp' = sp.current.shuff.means,
                              'sp.wtd' = sp.current.shuff.means))
                  
                } # shuffle.loop
      
      ## Close parallel backend
      stopCluster(cl)
      unregister_dopar()
      
      
      ## Reorganize
      # Separate shuff data by groups
      sp.shuff.data[[group.loop]] <- list()
      sp.shuff.data.wtd[[group.loop]] <- list()
      
      # Extract SP data
      for(sp.term.loop in terms){
        # sp.term.loop = terms[1]
        
        # Unweighted
        sp.shuff.data[[group.loop]][[sp.term.loop]] <- 
          lapply(shuff.means, function(x){
            return(x[['sp']][[sp.term.loop]])
          })
        # Combine into dataframe
        sp.shuff.data[[group.loop]][[sp.term.loop]] <-
          data.frame(t(bind_rows(sp.shuff.data[[group.loop]][[sp.term.loop]])))
        
        # Weighted
        if(! is.null(weights[[grouping.loop]])){ # If weights
          sp.shuff.data.wtd[[group.loop]][[sp.term.loop]] <- 
            lapply(shuff.means, function(x){
              return(x[['sp.wtd']][[sp.term.loop]])
            })
          # Combine into dataframe
          sp.shuff.data.wtd[[group.loop]][[sp.term.loop]] <-
            data.frame(t(bind_rows(sp.shuff.data.wtd[[group.loop]][[sp.term.loop]])))
        } # If weights
        
      }; rm(sp.term.loop)
      
    }; rm(group.loop)
    
    
    ### Get stats
    # Unweighted
    shuffled.assignments.mean <- lapply(sp.shuff.data, function(y){lapply(y, function(x){apply(x, 1, mean)})})
    shuffled.assignments.sd <- lapply(sp.shuff.data, function(y){lapply(y, function(x){apply(x, 1, sd)})})
    
    # Weighted
    if(! is.null(weights[[grouping.loop]])){ # If weights
      shuffled.assignments.mean.wtd <- lapply(sp.shuff.data.wtd, function(y){lapply(y, function(x){apply(x, 1, mean)})})
      shuffled.assignments.sd.wtd <- lapply(sp.shuff.data.wtd, function(y){lapply(y, function(x){apply(x, 1, sd)})})
    }else{
      shuffled.assignments.mean.wtd <- NA
      shuffled.assignments.sd.wtd <- NA
    } # If weights
    
    # Consecutive significant window minimum
    min.sig.window.ms <- 100
    
    ### Loop thru alphas
    sig.windows.shuffled.assignments <- list()
    sig.windows.shuffled.assignments.wtd <- list()
    for(alpha.loop in c(.05, .01, .001)){
      # alpha.loop = .05
      alpha.label <- paste0('alpha_',alpha.loop)
      
      ### Get CIs
      sp.ci <- list()
      
      ## One-tailed 
      alpha.z <- qnorm(1-alpha.loop) # one-tailed
      sp.ci <- list()
      sp.ci.wtd <- list()
      for(group.loop in names(groupings[[grouping.loop]])){
        # group.loop <- names(groupings[[grouping.loop]])[1]
        
        sp.ci[[group.loop]] <- list()
        sp.ci.wtd[[group.loop]] <- list()
        for(term.loop in terms){
          # term.loop <- terms[1]
          
          # Unweighted
          sp.ci[[group.loop]][[term.loop]] <- 
            shuffled.assignments.mean[[group.loop]][[term.loop]] + 
            (alpha.z * shuffled.assignments.sd[[group.loop]][[term.loop]])
          
          # Weighted
          if(! is.null(weights[[grouping.loop]])){ # If weights
            sp.ci.wtd[[group.loop]][[term.loop]] <- 
              shuffled.assignments.mean.wtd[[group.loop]][[term.loop]] + 
              (alpha.z * shuffled.assignments.sd.wtd[[group.loop]][[term.loop]])
          } # if weighted
        }; rm(term.loop)
      }; rm(group.loop)
      
      
      ### Get sig windows
      sig.windows.shuffled.assignments[[alpha.label]] <- list()
      sig.windows.shuffled.assignments.wtd[[alpha.label]] <- list()
      for(group.loop in names(groupings[[grouping.loop]])){
        # group.loop = names(groupings[[grouping.loop]])[1]
        
        sig.windows.shuffled.assignments[[alpha.label]][[group.loop]] <- list()
        sig.windows.shuffled.assignments.wtd[[alpha.label]][[group.loop]] <- list()
        
        ## Sig windows ("priming")
        for(term.loop in terms){
          # (term.loop = terms[1])
          
          # Unweighted
          current.data <- data.frame('mean' = sp.means[[group.loop]][[term.loop]],
                                     'upper.ci' = sp.ci[[group.loop]][[term.loop]])
          current.data$sig <- as.numeric((current.data$mean > current.data$upper.ci))
          sig.windows.shuffled.assignments[[alpha.label]][[group.loop]][[term.loop]] <- 
            get.significant.windows(current.data$sig,
                                    .sample.labels = row.names(current.data),
                                    include.duration = TRUE,
                                    .exclude.sig.durations.under.ms = min.sig.window.ms,
                                    .exclude.times.before.ms = -median.rt.times['sp'],
                                    output.class = "data.frame")
          rm(current.data)
          
          # Weighted
          if(! is.null(weights[[grouping.loop]])){ # If weights
            current.data <- data.frame('mean' = sp.means.wtd[[group.loop]][[term.loop]],
                                       'upper.ci' = sp.ci.wtd[[group.loop]][[term.loop]])
            current.data$sig <- as.numeric((current.data$mean > current.data$upper.ci))
            sig.windows.shuffled.assignments.wtd[[alpha.label]][[group.loop]][[term.loop]] <- 
              get.significant.windows(current.data$sig,
                                      .sample.labels = row.names(current.data),
                                      include.duration = TRUE,
                                      .exclude.sig.durations.under.ms = min.sig.window.ms,
                                      .exclude.times.before.ms = -median.rt.times['sp'],
                                      output.class = "data.frame")
            rm(current.data)
          } # if weights
          
        }; rm(term.loop)
        
      }; rm(group.loop)
      
      # Clean up
      rm(sp.ci, sp.ci.wtd)
      
    } # alpha.loop
    
    
    save.these <- c(
      'sig.windows.shuffled.assignments',
      'sig.windows.shuffled.assignments.wtd',
      'shuffled.assignments.mean',
      'shuffled.assignments.mean.wtd',
      'shuffled.assignments.sd',
      'shuffled.assignments.sd.wtd'
    )
    save.dir <- paste0(output.path,
                       '/data/1i - stats - shuffled cluster assignments/',
                       include.depth.dir,'/',
                       include.inactive.dir,'/',
                       "with hga/")
    dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
    save(list = save.these,
         file = paste0(save.dir,
                       grouping.loop,
                       ' - stats - shuffled cluster assignments.RData'))
    
    # Clean up
    rm(list = save.these)
    
  } # grouping.loop
  
# } # band.loop

# Finish!
message('Script completed successfully. ',Sys.time())







