### Get cluster significance by comparing clusters of RSA results on real vs. shuffled data
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
library('caret')
library('doParallel')
library('NMF') # for NMF
library('pals')
library('Hmisc')
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
    n.cores.to.use = 10
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 44
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
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))

## Function for weighted average
weighted.average <- function(x, w){
  ## Sum of the weights 
  sum.w <- sum(w, na.rm = TRUE)
  ## Sum of the weighted $x_i$ 
  xw <- sum(w*x, na.rm = TRUE)
  
  ## Return the weighted average 
  return(xw/sum.w)
}

## Function for getting weighted squared error (thx to Alex Stephenson)
weighted.se.mean <- function(x, w, na.rm = TRUE){
  if(na.rm){
    keep.these <- which(! is.na(x))
    x <- x[keep.these]
    w <- w[keep.these]
  }
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  numerator = sum(w*(x-weighted.average(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)
}


### Clean up
keep <- c(ls(), 
          'keep.all.this', 
          'band.loop', 
          'alpha.hps',
          'sig.hp.loop')

### Loop thru beta/high gamma data
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
      
      
      ##
      ## Load shuffle data
      ##
      
      # Paths
      shuffle.zs.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/shuffle zs and sig elecs for each shuffle loop/')
      shuffle.nmf.path <- paste0(output.path, 'data/1f - shuffled data nmf results/',
                                 include.depth.dir,'/',
                                 include.inactive.dir,'/',
                                 "with hga/",
                                 '/',sig.hp.loop,
                                 '/model weights/')
      # Files
      shuffle.zs.files <- list.files(shuffle.zs.path)
      shuffle.nmf.files <- list.files(shuffle.nmf.path)
      # Files to loop thru
      shuffle.loops.to.loop.thru <- 
        sort(as.numeric(intersect(gsub("zs and sig elecs - shuffle loop ","",
                                       gsub(".RData","", shuffle.zs.files, fixed = TRUE)),
                                  gsub("nmf - shuffled data - loop ","", 
                                       gsub(".RData","", shuffle.nmf.files, fixed = TRUE)))))
      
      # Read data
      ## Set up parallel processing
      # Close any old parallel backends
      unregister_dopar()
      # Set up parallel workers in case caret wants to use it
      cl <- makeCluster(n.cores.to.use, type = "FORK")
      registerDoParallel(cl)
      
      shuffle.data <- 
        foreach(shuffle.loop = shuffle.loops.to.loop.thru) %dopar% {
          # shuffle.loop = shuffle.loops.to.loop.thru[1]
          
          load(paste0(shuffle.zs.path, 'zs and sig elecs - shuffle loop ',shuffle.loop,'.RData')) # loads shuffle.zs.adjusted and shuffle.zs.sig.elecs
          load(paste0(shuffle.nmf.path, 'nmf - shuffled data - loop ',shuffle.loop,'.RData')) # loads shuffle.weights
          
          # Output
          return(list('shuffle.zs.adjusted' = shuffle.zs.adjusted,
                      'shuffle.weights' = shuffle.weights))
        } # shuffle.loop
      
      ## Close parallel backend
      stopCluster(cl)
      unregister_dopar()
      
      message("Done loading ",length(shuffle.data)," shuffle stats files! ",Sys.time())
      beep()
      
      ## Reorganize
      shuffle.weights <- lapply(shuffle.data, function(x){x[['shuffle.weights']]})
      shuffle.zs.adjusted <- lapply(shuffle.data, function(x){x[['shuffle.zs.adjusted']]})
      rm(shuffle.data)
      gc()
      
      
      ### Get shuffle.elec.weights
      shuffle.elec.weights <- list()
      for(shuffle.loop in 1:length(shuffle.weights)){
        # shuffle.loop = 1
        
        # Accidentally deleted "NMF_" in previous step... add back in
        if(! all(grepl("NMF_",shuffle.weights[[shuffle.loop]]$cluster))){
          shuffle.weights[[shuffle.loop]]$cluster <- paste0('NMF_',shuffle.weights[[shuffle.loop]]$cluster)
        } # If "NMF_" prefix missing
        
        shuffle.elec.weights[[shuffle.loop]] <- split(shuffle.weights[[shuffle.loop]], shuffle.weights[[shuffle.loop]]$cluster)
        
        for(cluster.loop in names(shuffle.elec.weights[[shuffle.loop]])){
          # cluster.loop = names(shuffle.elec.weights[[shuffle.loop]])[1]
          shuffle.elec.weights[[shuffle.loop]][[cluster.loop]]$weight <- 
            shuffle.elec.weights[[shuffle.loop]][[cluster.loop]][,paste0('weight_',cluster.loop)]
          shuffle.elec.weights[[shuffle.loop]][[cluster.loop]] <- unlist(t(shuffle.elec.weights[[shuffle.loop]][[cluster.loop]][,'weight', drop = FALSE]))
        }; rm(cluster.loop)
      }; rm(shuffle.loop)
      
      
      ##
      ## Load elec means
      ##
      
      message("Attach elec means!")
      attach(paste0(path,'/analysis/R/electrode means/means - warped data - multiband/output - ',band.loop,'/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses
      elec.means <- elec.means
      median.rt.samples <- median.rt.samples
      detach()
      message("Detach elec means!")
      
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
      
      
      ##
      ## Load adjusted t-values
      ##
      
      adjusted.zs.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
      adjusted.zs.file <- 'adjusted zs with sig windows and elecs.RData'
      load(paste0(adjusted.zs.path, adjusted.zs.file)) # loads significance.hps, zs.adjusted, zs.sig.elecs, and zs.sig.windows
      
      # Sample labels to plot for each task
      keep.sample.labels <- list()
      for(task.loop in names(median.rt.samples)){
        last.sample <- time.convert(500, "times", "samples")
        keep.sample.labels[[task.loop]] <- time.convert((-median.rt.samples[task.loop] - 100):last.sample,
                                                        "samples", "sample.labels")
        rm(last.sample)
      }; rm(task.loop)
      
      terms <- colnames(zs.adjusted[[1]])
      
      
      ## 
      ## Load cluster data
      ##
      
      rank.loop <- 5
      current.read.nmf.path <- paste0(output.path,'data/1e - nmf/',
                                      include.depth.dir,'/',
                                      include.inactive.dir,'/',
                                      "with hga/",
                                      sig.hp.loop,'/',
                                      'individual models/RData files/rank=',rank.loop,'/')
      message("Attaching NMF...")
      attach(paste0(current.read.nmf.path, 'NMF clustering results.RData')) # Loads: 'factors','weights','cluster.assignments','cluster.assignments.1hot','cluster.colors','cluster.colors.rgb','weights.and.clusters.localizations'
      weights <- weights
      cluster.assignments <- cluster.assignments
      detach()
      message("...detached!")
      
      # Get clusters, sorted
      clusters <- unique(cluster.assignments$cluster)
      clusters <- clusters[order(as.numeric(gsub("NMF_","",clusters)))]
      
      # Number of clusters
      n.clusters <- length(clusters) # might be different from rank.loop -- added 1 for subthreshhold elecs when trying out threshholding elecs by NMF weight for cluster membership
      
      # Get list of elecs per cluster
      cluster.elecs <- list()
      for(cluster.loop in clusters){
        cluster.elecs[[cluster.loop]] <- cluster.assignments[cluster.assignments$cluster == cluster.loop, 'electrode']
      }; rm(cluster.loop)
      
      
      ### Get cluster weights
      colnames(weights) <- paste0('weight_', colnames(weights))
      weights <- cbind(cluster.assignments[,c('electrode','cluster')], weights)
      elec.weights <- split(weights, weights$cluster)
      for(cluster.loop in clusters){
        # cluster.loop = clusters[1]
        elec.weights[[cluster.loop]]$weight <- elec.weights[[cluster.loop]][,paste0('weight_',cluster.loop)]
        elec.weights[[cluster.loop]] <- unlist(t(elec.weights[[cluster.loop]][,'weight',drop = FALSE]))
      }; rm(cluster.loop)
      
      
      ###
      ### Get cluster stats
      ###
      
      # Real data
      cluster.means <- list()
      wtd.cluster.means <- list()
      cluster.ses <- list()
      wtd.cluster.ses <- list()
      
      # Shuffle data
      shuffle.cluster.means <- list()
      wtd.shuffle.cluster.means <- list()
      shuffle.cluster.sds <- list()
      wtd.shuffle.cluster.sds <- list()
      
      for(cluster.loop in clusters){
        # cluster.loop = clusters[1]
        
        # Storage
        cluster.means[[cluster.loop]] <- list()
        cluster.ses[[cluster.loop]] <- list()
        shuffle.cluster.means[[cluster.loop]] <- list()
        shuffle.cluster.sds[[cluster.loop]] <- list()
        wtd.cluster.means[[cluster.loop]] <- list()
        wtd.cluster.ses[[cluster.loop]] <- list()
        wtd.shuffle.cluster.means[[cluster.loop]] <- list()
        wtd.shuffle.cluster.sds[[cluster.loop]] <- list()
        
        
        ##
        ## HGA
        ##
        
        for(task.loop in names(elec.means)){
          # (task.loop = names(elec.means)[3])
          
          # Get means and SEs
          cluster.means[[cluster.loop]][[task.loop]] <-
            apply(elec.means[[task.loop]][,cluster.elecs[[cluster.loop]]], 1, mean, na.rm = TRUE)
          cluster.ses[[cluster.loop]][[task.loop]] <-
            apply(elec.means[[task.loop]][,cluster.elecs[[cluster.loop]]], 1, function(x){
              x <- x[!is.na(x)] # manually get rid of NAs so SE works
              x <- sd(x) / sqrt(length(x))
              return(x)
            })
          
          # Weighted means and SEs
          wtd.cluster.means[[cluster.loop]][[task.loop]] <-
            apply(elec.means[[task.loop]][,cluster.elecs[[cluster.loop]]], 1, function(x){
              x <- weighted.average(x, w = elec.weights[[cluster.loop]][,cluster.elecs[[cluster.loop]]])
              return(x)
            })
          wtd.cluster.ses[[cluster.loop]][[task.loop]] <-
            apply(elec.means[[task.loop]][,cluster.elecs[[cluster.loop]]], 1, function(x){
              x <- weighted.se.mean(x, w = elec.weights[[cluster.loop]][,cluster.elecs[[cluster.loop]]])  
              return(x)
            })
          
        }; rm(task.loop)
        
        
        ## 
        ## Model terms
        ##
        
        for(term.loop in terms){
          # term.loop = terms[1]
          
          ### Real data
          current.data <- 
            do.call(cbind, lapply(zs.adjusted[cluster.elecs[[cluster.loop]]], function(x){
              x <- x[, term.loop, drop = FALSE]
              return(x)
            }))
          
          # Means and SEs
          cluster.means[[cluster.loop]][[term.loop]] <-
            apply(current.data, 1, mean)
          cluster.ses[[cluster.loop]][[term.loop]] <-
            apply(current.data, 1, sd) / sqrt(ncol(current.data))
          
          # Weighted means and SEs
          wtd.cluster.means[[cluster.loop]][[term.loop]] <-
            apply(current.data, 1, function(x){
              x <- weighted.average(x, w = elec.weights[[cluster.loop]][,cluster.elecs[[cluster.loop]]])
              return(x)
            })
          wtd.cluster.ses[[cluster.loop]][[term.loop]] <-
            apply(current.data, 1, function(x){
              x <- weighted.se.mean(x, w = elec.weights[[cluster.loop]][,cluster.elecs[[cluster.loop]]])
              return(x)
            })
          
          # Clean up
          rm(current.data)
          
          
          ### Shuffled data
          
          ## Get the shuffle cluster with the highest mean for this term
          current.means <- list()
          wtd.current.means <- list()
          
          for(shuffle.loop in 1:length(shuffle.zs.adjusted)){
            # shuffle.loop = 1
            
            current.means[[shuffle.loop]] <- list()
            wtd.current.means[[shuffle.loop]] <- list()
            
            for(shuffle.cluster.loop in clusters){
              # shuffle.cluster.loop = clusters[2]
              
              ### Cluster assignments based off shuffle-specific NMF clustering
              current.elecs <- colnames(shuffle.elec.weights[[shuffle.loop]][[shuffle.cluster.loop]])
              current.shuffle.data <- 
                do.call(cbind, lapply(shuffle.zs.adjusted[[shuffle.loop]][current.elecs], function(x){
                  x <- x[, term.loop, drop = FALSE]
                  return(x)
                }))
              # Store
              current.means[[shuffle.loop]][[shuffle.cluster.loop]] <- apply(current.shuffle.data, 1, mean)
              wtd.current.means[[shuffle.loop]][[shuffle.cluster.loop]] <- apply(current.shuffle.data, 1, function(x){
                x <- weighted.average(x, w = shuffle.elec.weights[[shuffle.loop]][[shuffle.cluster.loop]])
                return(x)
              })
              rm(current.shuffle.data)
              
            }; rm(shuffle.cluster.loop)
            
            # Combine
            current.means[[shuffle.loop]] <- do.call(cbind, current.means[[shuffle.loop]])
            wtd.current.means[[shuffle.loop]] <- do.call(cbind, wtd.current.means[[shuffle.loop]])
            
            # Get the cluster with the max for this term
            current.means[[shuffle.loop]] <- 
              current.means[[shuffle.loop]][,names(which.max(colMeans(current.means[[shuffle.loop]])))[1]]
            wtd.current.means[[shuffle.loop]] <- 
              wtd.current.means[[shuffle.loop]][,names(which.max(colMeans(wtd.current.means[[shuffle.loop]])))[1]]
            
          }; rm(shuffle.loop)
          
          ## Shuffle-loop specific NMF
          current.means <- do.call(cbind, current.means)
          wtd.current.means <- do.call(cbind, wtd.current.means)
          shuffle.cluster.means[[cluster.loop]][[term.loop]] <- 
            apply(current.means, 1, mean)
          shuffle.cluster.sds[[cluster.loop]][[term.loop]] <-
            apply(current.means, 1, sd)
          wtd.shuffle.cluster.means[[cluster.loop]][[term.loop]] <- 
            apply(wtd.current.means, 1, mean)
          wtd.shuffle.cluster.sds[[cluster.loop]][[term.loop]] <- 
            apply(wtd.current.means, 1, sd)
          
          # Clean up
          rm(current.means)
          
        }; rm(term.loop)
        
        print(paste0("Summary stats done for cluster loop ",cluster.loop,". ", Sys.time()))
      }; rm(cluster.loop)
      
      
      ###
      ### Get sig windows
      ###
      
      # Sig window storage
      sig.cluster.windows <- list()
      wtd.sig.cluster.windows <- list()
      
      # Confidence interval storage
      cis <- list()
      wtd.cis <- list()
      
      for(alpha.loop in c(.05, .01, .001)){
        # alpha.loop = .01
        
        # Get metadata
        alpha.label <- paste0("alpha=",alpha.loop)
        alpha.z <- qnorm(1 - alpha.loop) # 1-tailed
        
        # Sig window storage
        sig.cluster.windows[[alpha.label]] <- list()
        wtd.sig.cluster.windows[[alpha.label]] <- list()
        
        # Confidence interval storage
        cis[[alpha.label]] <- list()
        wtd.cis[[alpha.label]] <- list()
        
        ## Loop thru clusters
        for(cluster.loop in clusters){
          # cluster.loop = clusters[1]
          
          # Sig window storage
          sig.cluster.windows[[alpha.label]][[cluster.loop]] <- list()
          wtd.sig.cluster.windows[[alpha.label]][[cluster.loop]] <- list()
          
          # Confidence interval storage
          cis[[alpha.label]][[cluster.loop]] <- list()
          wtd.cis[[alpha.label]][[cluster.loop]] <- list()
          
          ## Loop thru terms
          for(term.loop in terms){
            # term.loop = terms[1]
            
            # CIs
            cis[[alpha.label]][[cluster.loop]] <- 
              (alpha.z * shuffle.cluster.sds[[cluster.loop]][[term.loop]])
            wtd.cis[[alpha.label]][[cluster.loop]] <- 
              (alpha.z * wtd.shuffle.cluster.sds[[cluster.loop]][[term.loop]])
            
            # Each sample sig?
            sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]] <-
              as.numeric(cluster.means[[cluster.loop]][[term.loop]] > 
                           (shuffle.cluster.means[[cluster.loop]][[term.loop]] + cis[[alpha.label]][[cluster.loop]]))
            wtd.sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]] <-
              as.numeric(wtd.cluster.means[[cluster.loop]][[term.loop]] > 
                           (wtd.shuffle.cluster.means[[cluster.loop]][[term.loop]] + wtd.cis[[alpha.label]][[cluster.loop]]))
            
            # Get sig windows
            sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]] <-
              get.significant.windows(sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]],
                                      .sample.labels = names(shuffle.cluster.means[[cluster.loop]][[term.loop]]),
                                      .exclude.times.before.ms = -time.convert(median.rt.samples['sp'], "samples", "times"),
                                      .exclude.sig.durations.under.ms = 100,
                                      output.class = 'data.frame',
                                      include.duration = TRUE)
            wtd.sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]] <-
              get.significant.windows(wtd.sig.cluster.windows[[alpha.label]][[cluster.loop]][[term.loop]],
                                      .sample.labels = names(wtd.shuffle.cluster.means[[cluster.loop]][[term.loop]]),
                                      .exclude.times.before.ms = -time.convert(median.rt.samples['sp'], "samples", "times"),
                                      .exclude.sig.durations.under.ms = 100,
                                      output.class = 'data.frame',
                                      include.duration = TRUE)
            
          }; rm(term.loop)
        }; rm(cluster.loop)
        
        rm(alpha.label, alpha.z)
      }; rm(alpha.loop)
      
      
      
      ### Reorganize
      cluster.means <- list(
        'un' = cluster.means,
        'wtd' = wtd.cluster.means
      )
      cluster.ses <- list(
        'un' = cluster.ses,
        'wtd' = wtd.cluster.ses
      )
      shuffle.cluster.means <- list(
        'un' = shuffle.cluster.means,
        'wtd' = wtd.shuffle.cluster.means
      )
      shuffle.cluster.sds <- list(
        'un' = shuffle.cluster.sds,
        'wtd' = wtd.shuffle.cluster.sds
      )
      shuffle.cluster.sig.windows <- list(
        'un' = sig.cluster.windows,
        'wtd' = wtd.sig.cluster.windows
      )
      shuffle.cluster.cis <- list(
        'un' = cis,
        'wtd' = wtd.cis
      )
      
      
      ### Save everything
      
      ###
      ### Save
      ###
      
      save.these <- c(
        # 'shuffle.median',
        # 'shuffle.025th.quantile',
        # 'shuffle.95th.quantile',
        # 'shuffle.975th.quantile',
        # 'shuffle.99th.quantile',
        # 'shuffle.999th.quantile',
        'cluster.means',
        'cluster.ses',
        'cluster.elecs',
        'shuffle.cluster.means',
        'shuffle.cluster.sds',
        'shuffle.cluster.sig.windows',
        'shuffle.cluster.cis'
      )
      
      save.stats.path <- paste0(output.path, 'data/1g - cluster stats - cluster means vs means from clusters of shuffled data/',
                                include.depth.dir,'/',
                                include.inactive.dir,'/',
                                "with hga/",
                                '/',sig.hp.loop,'/')
      dir.create(save.stats.path, showWarnings = FALSE, recursive = TRUE)
      save(list = save.these,
           file = paste0(save.stats.path, 'all stats.RData'))
      
    # } # sig.hp.loop
    
    
    # } # band.loop
    
    # Finish!
    message('Script completed successfully. ',Sys.time())
    
    
    
    
    
    
    
    