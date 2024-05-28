### Shuffle elecs across clusters to get noise distribution (conservative since it's only significant elecs)
### adam.milton.morgan@gmail.com

###
### Readme
###

## Tried:
# - shuff analysis
# - chi squared across clusters (doesn't capture MFG)
# - chi squared across regions (doesn't capture IFG)
# - define how many you'd expect in each region given proportion significant overall evenly distributed across regions (excluding current region so as not to bias), then simulate (sample 1s -- in cluster -- and 0s -- not in cluster with corresponding probabilities (proportion elecs in cluster * proportion elecs significant), and compare actual to distribution. one-tailed tests without corrections give you IFG and MFG in late syntax (neither survives corrections) and MFG in early syntax (survives).  So... not ideal.

## To try:
# I think the way forward might be two kinds of chi squared tests: (1) across clusters, to ask whether IFG is above chance for late cluster, and (2) across regions, to ask whether MFG (and maybe combine STG and MTG?) are more represented than other areas in each cluster. Will have to correct for multiple comparisons, but I think these will survive


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
    n.cores.to.use = 12
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
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Violin plots
source(paste0(path,'/analysis/R/functions/plot_violin.R'))
# Line chart
source(paste0(path,'/analysis/R/functions/plot_line_chart.R'))


# Colors for plotting
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

# Regions
gross.rois <- read.csv(paste0(path,'analysis/R/define ROIs/output/gross ROIs.csv'))
rownames(gross.rois) <- gross.rois$region.clinical
gross.rois[gross.rois$region.category == 'pericentral','region.category'] <- 'SMC'
gross.rois[gross.rois$grosser == 'pericentral','grosser'] <- 'SMC'

### Clean up
keep.all.this <- c(ls(), 
                   'keep.all.this', 
                   'band.loop')

### Loop thru beta/high gamma data
# for(band.loop in c('high_gamma','beta')){
band.loop = c('high_gamma','beta')[1]

      rm(list = ls()[! ls() %in% keep.all.this])
      gc()
      
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
      elec.info$gross.roi <- gross.rois[elec.info$region_clinical,'region.category']
      
      
        ## Save directory
      output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/')
      
      
      ### Load electrode importances
      elec.zs.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
      elec.zs.file <- 'adjusted zs with sig windows and elecs.RData'
      load(paste0(elec.zs.path, elec.zs.file))
      
      ### RTs
      load(paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
      median.rt.times <- time.convert(median.rt.samples, "samples", "times")
      
      
      ### Load elec means for plotting here and below
      ## Load non-band.loop means
      message('Attaching elec means from non-current band.loop...')
      attach(paste0(path,'analysis/R/electrode means/means - warped data - multiband/output - ',
                    ifelse(band.loop == 'beta','high_gamma','beta'),
                    '/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses
      other.band.elec.means <- elec.means
      detach()
      message('...detached!')
      # Add NA dataframe for Patient006 LP block
      blank.Patient006.lp.df <- data.frame(matrix(nrow = nrow(other.band.elec.means$lp[[1]]),
                                             ncol = ncol(other.band.elec.means$pn[['Patient006']])))
      colnames(blank.Patient006.lp.df) <- colnames(other.band.elec.means$pn[['Patient006']])
      rownames(blank.Patient006.lp.df) <- rownames(other.band.elec.means$lp[[1]])
      other.band.elec.means$lp[['Patient006']] <- blank.Patient006.lp.df
      
      # Collapse across patients
      other.band.elec.means <- lapply(other.band.elec.means, bind_cols)
      
      
      ## Load band.loop means
      load(paste0(path,'analysis/R/electrode means/means - warped data - multiband/output - ',band.loop,'/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses
      rm(elec.sds, elec.ses)
      
      # Get metadata
      tasks <- names(elec.means)
      
      # Add NA dataframe for Patient006 LP block
      blank.Patient006.lp.df <- data.frame(matrix(nrow = nrow(elec.means$lp[[1]]),
                                             ncol = ncol(elec.means$pn[['Patient006']])))
      colnames(blank.Patient006.lp.df) <- colnames(elec.means$pn[['Patient006']])
      rownames(blank.Patient006.lp.df) <- rownames(elec.means$lp[[1]])
      elec.means$lp[['Patient006']] <- blank.Patient006.lp.df
      
      # Collapse across patients
      elec.means <- lapply(elec.means, bind_cols)
      
      
        # Significance hyperparameters to loop through
        alpha.hps <- c('alpha=.05_min=100ms',
                       'alpha=.01_min=50ms')
        
        # for(alpha.loop in alpha.hps){ ## UNCOMMENT
          alpha.loop = alpha.hps[1]
          
          print(paste0("Beginning loop: ",
                       alpha.loop,"."))
          
          ## Path info
          nmf.results.path <- paste0(output.path,'data/1t - clustering just syntax/',
                                     alpha.loop,'/',
                                     'all models/')
          nmf.results.file <- 'NMF models of syntax elecs.RData'
          
          ## Define ROIs
          rois.to.plot <- c('IFG','MFG','SMC','STG','MTG','IPL')
          
          ## Select appropriate elecs depending on loops
          all.roi.elecs.analyzed <- elec.info[(elec.info$gross.roi %in% rois.to.plot) & 
                                                (elec.info$use.these.elecs == 1),
                                              c('patient_elec','gross.roi')]
          elecs.with.syntax <- zs.sig.elecs[[alpha.loop]]$diff.voice
          all.roi.elecs.analyzed$syntax.sig <- 
            as.numeric(all.roi.elecs.analyzed$patient_elec %in% elecs.with.syntax)
          
          ### Load models
          load(paste0(nmf.results.path, nmf.results.file))
          
          # Model ranks
          nmf.ranks.to.try <- nmf.models$measures$rank
          
          
          ### Cluster data for line plots and brain plots
          ## Loop thru NMFS and create data for brain plots
          for(nmf.rank.loop in nmf.ranks.to.try){ # UNCOMMENT
            # nmf.rank.loop = nmf.ranks.to.try[2]
            current.label = paste0('rank=',nmf.rank.loop)  
            
            # Get current cluster assignments dataframe
            cluster.assignments <-
              read.csv(paste0(output.path,'data/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/csvs for brain plots/rank=',nmf.rank.loop,'/cluster_assignments.csv'))
            row.names(cluster.assignments) <- cluster.assignments$electrode
            cluster.assignments$gross.roi <- elec.info[cluster.assignments$electrode, 'gross.roi']
            
            # Get NMF weights
            cluster.assignments <-
              cbind(cluster.assignments,
                    read.csv(paste0(output.path,'data/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/csvs for brain plots/rank=',nmf.rank.loop,'/weights.csv')))
            
            # Add column for winning weight
            cluster.assignments$winning.weight <- NA
            for(elec.loop in 1:nrow(cluster.assignments)){
              cluster.assignments[elec.loop, 'winning.weight'] <- cluster.assignments[elec.loop, cluster.assignments$cluster[elec.loop]]
            }; rm(elec.loop)
            
            # Clusters
            clusters <- sort(unique(cluster.assignments$cluster))
            
            ### Get means and standard errors
            ## Store values for plots
            cluster.elecs <- list()
            cluster.elec.means.weighted <- list()
            cluster.elec.ses.weighted <- list()
            cluster.elec.maxes.weighted <- list()
            cluster.elec.means.unweighted <- list()
            cluster.elec.ses.unweighted <- list()
            cluster.elec.maxes.unweighted <- list()
            cluster.xs <- list()
            cluster.z.means.weighted <- list()
            cluster.z.ses.weighted <- list()
            cluster.z.maxes.weighted <- list()
            cluster.z.means.unweighted <- list()
            cluster.z.ses.unweighted <- list()
            cluster.z.maxes.unweighted <- list()
            cluster.z.xs <- list()
            for(cluster.loop in clusters){
              # cluster.loop = clusters[1]
              
              ### Define functions with current weights
              # Current weights
              cluster.weights <-
                lapply(split(cluster.assignments, cluster.assignments$cluster), 
                       function(x){return(x$winning.weight)})[[cluster.loop]]
              
              ## Function for weighted average with these weights (cluster-specific!)
              weighted.average <- function(x){
                ## Remove NAs
                keep.these <- which(! is.na(x))
                x <- x[keep.these]
                fn.cluster.weights <- cluster.weights[keep.these]
                
                ## Sum of the weights 
                sum.w <- sum(fn.cluster.weights)
                ## Sum of the weighted $x_i$ 
                xw <- sum(fn.cluster.weights*x)
                
                ## Return the weighted average 
                return(xw/sum.w)
              } # weighted.average()
              
              ## Function for getting weighted squared error (thx to Alex Stephenson) with these weights (cluster-specific!)
              weighted.se.mean <- function(x){
                ## Remove NAs
                keep.these <- which(! is.na(x))
                x <- x[keep.these]
                fn.cluster.weights <- cluster.weights[keep.these]
                
                
                ## Calculate effective N and correction factor
                n_eff <- (sum(fn.cluster.weights))^2/(sum(fn.cluster.weights^2))
                correction = n_eff/(n_eff-1)
                
                ## Get weighted variance 
                numerator = sum(fn.cluster.weights*(x-weighted.average(x))^2)
                denominator = sum(fn.cluster.weights)
                
                ## get weighted standard error of the mean 
                se_x = sqrt((correction * (numerator/denominator))/n_eff)
                return(se_x)
              } # weighted.se.mean()
              
              ## Function for getting unweighted standard error
              se.mean <- function(x, na.rm = TRUE){
                if(na.rm){
                  keep.these <- which(! is.na(x))
                  x <- x[keep.these]
                }
                se <- sd(x) / sqrt(length(x))
                return(se)
              } # se.mean()
              
              
              ### Get summary stats
              # Get elecs in this cluster
              cluster.elecs[[cluster.loop]] <-
                lapply(split(cluster.assignments, cluster.assignments$cluster), function(x){return(x$electrode)})[[cluster.loop]]
              
              # ECoG vals, weighted
              cluster.elec.means.weighted[[cluster.loop]] <-
                lapply(elec.means, function(x){apply(x[,cluster.elecs[[cluster.loop]]], 1, weighted.average)})
              cluster.elec.ses.weighted[[cluster.loop]] <-
                lapply(elec.means, function(x){apply(x[,cluster.elecs[[cluster.loop]]], 1, weighted.se.mean)})
              cluster.elec.maxes.weighted[[cluster.loop]] <-
                max(unlist(cluster.elec.means.weighted[[cluster.loop]]) + unlist(cluster.elec.ses.weighted[[cluster.loop]]))
              cluster.xs[[cluster.loop]] <- lapply(cluster.elec.means.weighted[[cluster.loop]], names)
              
              # ECoG vals, unweighted
              cluster.elec.means.unweighted[[cluster.loop]] <-
                lapply(elec.means, function(x){apply(x[,cluster.elecs[[cluster.loop]]], 1, mean, na.rm = TRUE)})
              cluster.elec.ses.unweighted[[cluster.loop]] <-
                lapply(elec.means, function(x){apply(x[,cluster.elecs[[cluster.loop]]], 1, se.mean)})
              cluster.elec.maxes.unweighted[[cluster.loop]] <-
                max(unlist(cluster.elec.means.unweighted[[cluster.loop]]) + unlist(cluster.elec.ses.unweighted[[cluster.loop]]))
              
              # RSA vals, weighted
              cluster.z.means.weighted[[cluster.loop]] <-
                as.list(data.frame(elementwise.matrix.apply(zs.adjusted[cluster.elecs[[cluster.loop]]], .function = "weighted.average")))
              cluster.z.ses.weighted[[cluster.loop]] <-
                as.list(data.frame(elementwise.matrix.apply(zs.adjusted[cluster.elecs[[cluster.loop]]], .function = "weighted.se.mean")))
              cluster.z.maxes.weighted[[cluster.loop]] <-
                max(unlist(cluster.z.means.weighted[[cluster.loop]]) + unlist(cluster.z.ses.weighted[[cluster.loop]]))
              cluster.z.xs[[cluster.loop]] <- rownames(zs.adjusted[[1]])
              
              # RSA vals, unweighted
              cluster.z.means.unweighted[[cluster.loop]] <-
                as.list(data.frame(elementwise.matrix.apply(zs.adjusted[cluster.elecs[[cluster.loop]]], .function = "mean")))
              cluster.z.ses.unweighted[[cluster.loop]] <-
                as.list(data.frame(elementwise.matrix.apply(zs.adjusted[cluster.elecs[[cluster.loop]]], .function = "se.mean")))
              cluster.z.maxes.unweighted[[cluster.loop]] <-
                max(unlist(cluster.z.means.unweighted[[cluster.loop]]) + unlist(cluster.z.ses.unweighted[[cluster.loop]]))
              
              ## Set HGA at samples before stimulus to NA
              hga.sample.labels <- lapply(cluster.elec.means.weighted[[cluster.loop]], names)
              hga.times <- lapply(hga.sample.labels, time.convert, "sample.labels", "times")
              ecog.times <- time.convert(cluster.z.xs[[cluster.loop]], "sample.labels", "times")
              
              for(task.loop in names(cluster.elec.means.weighted[[cluster.loop]])){
                # task.loop = names(cluster.elec.means.weighted[[cluster.loop]])[2]
                
                # NA indices
                cluster.ecog.na.indices <- 
                  which(hga.times[[task.loop]] < -median.rt.times[task.loop] - 150)
                
                # Remove pre-stim HGA data, weighted
                cluster.elec.means.weighted[[cluster.loop]][[task.loop]] <-
                  cluster.elec.means.weighted[[cluster.loop]][[task.loop]][-cluster.ecog.na.indices]
                cluster.elec.ses.weighted[[cluster.loop]][[task.loop]] <-
                  cluster.elec.ses.weighted[[cluster.loop]][[task.loop]][-cluster.ecog.na.indices]
                cluster.xs[[cluster.loop]][[task.loop]] <-
                  cluster.xs[[cluster.loop]][[task.loop]][-cluster.ecog.na.indices]
                
                # Remove pre-stim HGA data, unweighted
                cluster.elec.means.unweighted[[cluster.loop]][[task.loop]] <-
                  cluster.elec.means.unweighted[[cluster.loop]][[task.loop]][-cluster.ecog.na.indices]
                cluster.elec.ses.unweighted[[cluster.loop]][[task.loop]] <-
                  cluster.elec.ses.unweighted[[cluster.loop]][[task.loop]][-cluster.ecog.na.indices]
                
                # Remove pre-stim RSA data
                if(task.loop == 'sp'){
                  # NA indices
                  cluster.rsa.na.indices <- 
                    which(ecog.times < -median.rt.times[task.loop] - 150)
                  
                  # Remove pre-stim RSA data, unweighted
                  cluster.z.means.weighted[[cluster.loop]] <-
                    lapply(cluster.z.means.weighted[[cluster.loop]], function(x){return(x[-cluster.rsa.na.indices])})
                  cluster.z.ses.weighted[[cluster.loop]] <-
                    lapply(cluster.z.ses.weighted[[cluster.loop]], function(x){return(x[-cluster.rsa.na.indices])})
                  cluster.z.xs[[cluster.loop]] <-
                    cluster.z.xs[[cluster.loop]][-cluster.rsa.na.indices]
                  
                  # Remove pre-stim RSA data, unweighted
                  cluster.z.means.unweighted[[cluster.loop]] <-
                    lapply(cluster.z.means.unweighted[[cluster.loop]], function(x){return(x[-cluster.rsa.na.indices])})
                  cluster.z.ses.unweighted[[cluster.loop]] <-
                    lapply(cluster.z.ses.unweighted[[cluster.loop]], function(x){return(x[-cluster.rsa.na.indices])})
                  
                  # Clean up 
                  rm(cluster.rsa.na.indices)
                } # if(task.loop == 'sp')
                
                # Clean up
                rm(cluster.ecog.na.indices)
              }; rm(task.loop)
              
              
              # Clean up
              rm(cluster.weights, weighted.average, weighted.se.mean)
              
            }; rm(cluster.loop)
            
            
            ### Get shuffle distributions
            ## Do it without weights -- too confusing otherwise
            n.shuffles <- 10000
            
            set.seed(123)
            shuffle.elec.means.temp <- list()
            shuffle.z.means.temp <- list()
            for(shuffle.loop in 1:n.shuffles){
              # shuffle.loop = 1
              
              # Shuffle cluster assignments
              shuffled.assignments <- data.frame('electrode' = cluster.assignments$electrode,
                                            'cluster' = sample(cluster.assignments$cluster))
              shuffled.cluster.means <- list()
              
              # Get cluster means by task (ECoG) and term (RSA)
              shuffle.elec.means.temp[[shuffle.loop]] <- list()
              shuffle.z.means.temp[[shuffle.loop]] <- list()
              for(cluster.loop in clusters){
                # cluster.loop = clusters[1]
                
                # Get elecs in this cluster
                .current.elecs <-
                  lapply(split(shuffled.assignments, shuffled.assignments$cluster), function(x){return(x$electrode)})[[cluster.loop]]
                
                # ECoG vals
                shuffle.elec.means.temp[[shuffle.loop]][[cluster.loop]] <-
                  lapply(elec.means, function(x){apply(x[,.current.elecs], 1, mean, na.rm = TRUE)})
                
                # RSA vals, unweighted
                mean.with.na.rm <- function(x){mean(x, na.rm = TRUE)}
                shuffle.z.means.temp[[shuffle.loop]][[cluster.loop]] <-
                  as.list(data.frame(elementwise.matrix.apply(zs.adjusted[.current.elecs], .function = "mean.with.na.rm")))
                shuffle.z.means.temp[[shuffle.loop]][[cluster.loop]] <-
                  lapply(shuffle.z.means.temp[[shuffle.loop]][[cluster.loop]], function(x){
                    names(x) <- rownames(zs.adjusted[[1]])
                    return(x)
                  })
                
                ## Set HGA at samples before stimulus to NA
                keep.hga.sample.labels <- cluster.xs[[cluster.loop]]
                keep.rsa.sample.labels <- cluster.z.xs[[cluster.loop]]
                
                for(task.loop in names(shuffle.elec.means.temp[[shuffle.loop]][[cluster.loop]])){
                  # task.loop = names(shuffle.elec.means.temp[[shuffle.loop]][[cluster.loop]])[1]
                  
                  # Remove pre-stim HGA data, weighted
                  shuffle.elec.means.temp[[shuffle.loop]][[cluster.loop]][[task.loop]] <-
                    shuffle.elec.means.temp[[shuffle.loop]][[cluster.loop]][[task.loop]][keep.hga.sample.labels[[task.loop]]]
                  
                  # Remove pre-stim RSA data
                  if(task.loop == 'sp'){
                    # Remove pre-stim RSA data
                    shuffle.z.means.temp[[shuffle.loop]][[cluster.loop]] <-
                      lapply(shuffle.z.means.temp[[shuffle.loop]][[cluster.loop]], function(x){return(x[keep.rsa.sample.labels])})
                  } # if(task.loop == 'sp')
                  
                }; rm(task.loop)
                
              }#; rm(cluster.loop)
              
            }#; rm(shuffle.loop)
            
            
            ### Reorganize shuffled data
            shuffle.elec.means <- list()
            shuffle.z.means <- list()
            for(cluster.loop in clusters){
              # cluster.loop = clusters[1]
              
              shuffle.elec.means[[cluster.loop]] <- list()
              shuffle.z.means[[cluster.loop]] <- list()
              
              ### ECoG data (tasks)
              for(task.loop in tasks){
                # task.loop = tasks[1]
                shuffle.elec.means[[cluster.loop]][[task.loop]] <- 
                  data.frame(bind_rows(lapply(shuffle.elec.means.temp, function(x){x[[cluster.loop]][[task.loop]]})))
              }#; rm(task.loop)
              
              ### RSA data (terms)
              for(term.loop in names(cluster.z.means.unweighted[[cluster.loop]])){
                # term.loop = names(cluster.z.means.unweighted[[cluster.loop]])[1]
                shuffle.z.means[[cluster.loop]][[term.loop]] <-
                  data.frame(bind_rows(lapply(shuffle.z.means.temp, function(x){x[[cluster.loop]][[term.loop]]})))
              }#; rm(term.loop)
              
            }#; rm(cluster.loop)
            # rm(shuffle.elec.means.temp, shuffle.z.means.temp)
            
            
            ### Reorganize real data into dataframes with p-values etc.
            cluster.elec.stats <- list()
            cluster.elec.sig.windows <- list()
            cluster.z.stats <- list()
            cluster.z.sig.windows <- list()
            for(cluster.loop in clusters){
              # cluster.loop = clusters[1]
              
              ### ECoG
              cluster.elec.stats[[cluster.loop]] <- list()
              cluster.elec.sig.windows[[cluster.loop]] <- list()
              for(task.loop in tasks){
                # task.loop = tasks[1]
                
                cluster.elec.stats[[cluster.loop]][[task.loop]] <- data.frame(
                  'sample.label' = names(cluster.elec.means.weighted[[cluster.loop]][[task.loop]]),
                  'mean.wtd' = cluster.elec.means.weighted[[cluster.loop]][[task.loop]],
                  'se.wtd' = cluster.elec.ses.weighted[[cluster.loop]][[task.loop]],
                  'mean.un' = cluster.elec.means.unweighted[[cluster.loop]][[task.loop]],
                  'se.un' = cluster.elec.ses.unweighted[[cluster.loop]][[task.loop]])
                cluster.elec.stats[[cluster.loop]][[task.loop]]$p <- NA
                for(sample.loop in cluster.elec.stats[[cluster.loop]][[task.loop]]$sample.label){
                  # sample.loop = cluster.elec.stats[[cluster.loop]][[task.loop]]$sample.label[1]
                  cluster.elec.stats[[cluster.loop]][[task.loop]][sample.loop,'p'] <-
                    mean(as.numeric(cluster.elec.stats[[cluster.loop]][[task.loop]][sample.loop,'mean.un'] < 
                                      shuffle.elec.means[[cluster.loop]][[task.loop]][,sample.loop]))
                }; rm(sample.loop)
                  
                ## Sig windows
                # Default
                alpha = .05
                min.ms.sig.threshold = 100
                if(alpha.loop == "alpha=.01_min=50ms"){
                  alpha = .01
                  min.ms.sig.threshold = 50
                } # if(alpha.loop == "alpha=.01_min=50ms"){
                cluster.elec.sig.windows[[cluster.loop]][[task.loop]] <-
                  get.significant.windows(sig.vals = as.numeric(cluster.elec.stats[[cluster.loop]][[task.loop]]$p < alpha),
                                          .sample.labels = cluster.elec.stats[[cluster.loop]][[task.loop]]$sample.label,
                                          output.class = 'data.frame',include.duration = TRUE,
                                          .exclude.times.before.ms = -median.rt.times[task.loop],
                                          .exclude.sig.durations.under.ms = min.ms.sig.threshold)
                
              }#; rm(task.loop)
              
              ### RSA
              cluster.z.stats[[cluster.loop]] <- list()
              cluster.z.sig.windows[[cluster.loop]] <- list()
              for(term.loop in names(shuffle.z.means[[cluster.loop]])){
                # term.loop = names(shuffle.z.means[[cluster.loop]])[1]
                cluster.z.stats[[cluster.loop]][[term.loop]] <- data.frame(
                  'sample.label' = names(shuffle.z.means[[cluster.loop]][[term.loop]]),
                  'mean.wtd' = cluster.z.means.weighted[[cluster.loop]][[term.loop]],
                  'se.wtd' = cluster.z.ses.weighted[[cluster.loop]][[term.loop]],
                  'mean.un' = cluster.z.means.unweighted[[cluster.loop]][[term.loop]],
                  'se.un' = cluster.z.ses.unweighted[[cluster.loop]][[term.loop]])
                rownames(cluster.z.stats[[cluster.loop]][[term.loop]]) <-
                  cluster.z.stats[[cluster.loop]][[term.loop]]$sample.label
                cluster.z.stats[[cluster.loop]][[term.loop]]$p <- NA
                for(sample.loop in cluster.z.stats[[cluster.loop]][[term.loop]]$sample.label){
                  # sample.loop = cluster.z.stats[[cluster.loop]][[term.loop]]$sample.label[1]
                  cluster.z.stats[[cluster.loop]][[term.loop]][sample.loop,'p'] <-
                    mean(as.numeric(cluster.z.stats[[cluster.loop]][[term.loop]][sample.loop,'mean.un'] < 
                                      shuffle.z.means[[cluster.loop]][[term.loop]][,sample.loop]))
                }; rm(sample.loop)
                
                ## Sig windows
                cluster.z.sig.windows[[cluster.loop]][[term.loop]] <-
                  get.significant.windows(sig.vals = as.numeric(cluster.z.stats[[cluster.loop]][[term.loop]]$p < alpha),
                                          .sample.labels = cluster.z.stats[[cluster.loop]][[term.loop]]$sample.label,
                                          output.class = 'data.frame',include.duration = TRUE,
                                          .exclude.times.before.ms = -median.rt.times['sp'],
                                          .exclude.sig.durations.under.ms = min.ms.sig.threshold)
                
              }#; rm(term.loop)
              
            }#; rm(cluster.loop)
            
            cluster.elec.sig.windows    
            cluster.z.sig.windows
            
            save.these <- c(
              'cluster.elec.stats',
              'cluster.z.stats',
              'cluster.elec.sig.windows',    
              'cluster.z.sig.windows',
              'cluster.elec.maxes.weighted',
              'cluster.elec.maxes.unweighted',
              'cluster.z.maxes.weighted',
              'cluster.z.maxes.unweighted'
            )
            save.dir <- paste0(output.path, 'data/1t - clustering just syntax/1t_b - values and stats for plots/',
                                                    alpha.loop,'/',
                                                    'rank=',nmf.rank.loop,'/')
            dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
            save(list = save.these,
                 file = paste0(save.dir, 'syntax clusters - summary stats and plot values.RData'))
            
            ## Clean up
            rm(list = save.these)
            gc()
            
            }; rm(nmf.rank.loop)
        } # alpha.loop
# }; rm(band.loop)



# Finish!
message('Script completed successfully. ',Sys.time())









