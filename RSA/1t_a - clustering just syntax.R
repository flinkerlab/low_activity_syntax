### Temporal clustering of syntax electrodes
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


### Clean up
keep.all.this <- c(ls(), 
                   'keep.all.this', 
                   'band.loop')

### Loop thru beta/high gamma data
# for(band.loop in c('high_gamma','beta')){ ## UNCOMMENT
band.loop = c('high_gamma','beta')[1] ## UNCOMMENT


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


use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec  


## Save directory
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/')


### Load electrode importances
elec.zs.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
elec.zs.file <- 'adjusted zs with sig windows and elecs.RData'
load(paste0(elec.zs.path, elec.zs.file))


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


###
### NMF Clustering
###


# Significance hyperparameters to loop through
alpha.hps <- c('alpha=.05_min=100ms',
               'alpha=.01_min=50ms')

for(alpha.loop in alpha.hps){
  # alpha.loop = alpha.hps[1]
  
  print(paste0("Beginning loop: ",
               alpha.loop,"."))
  
  # NMF ranks to try
  nmf.ranks.to.try <- 2:7
  
  ## Path info
  nmf.results.path <- paste0(output.path,'data/1t - clustering just syntax/',
                             alpha.loop,'/',
                             'all models/')
  nmf.results.file <- 'NMF models of syntax elecs.RData'
  
  ### IF already run, just load
  if(file.exists(paste0(nmf.results.path, nmf.results.file))){
    message("NMF previously run! Rename or delete old output to run this scrip for updated data output. (Script will update plots and CSVs for brain plots regardless.)")
    load(paste0(nmf.results.path, nmf.results.file))
  }else{ # if(file.exists(paste0(nmf.results.path, nmf.results.file))){
    
    ##
    ## Stack data
    ##
    
    # Initialize storage
    stacked.data <- list()
    
    
    ### Only elecs with significant syntax
    # Get syntax sig windows for this alpha.loop
    sig.syntax.windows <- lapply(zs.sig.windows, function(x){return(x$diff.voice[[alpha.loop]])})
    # Remove elecs only sig after production onset
    lapply(sig.syntax.windows, function(x){
      if(nrow(x) > 0){
        # Set last sig time to 0
        x$end.time <- sapply(x$end.time, function(y){min(y, 0)})
        # Update sig durations
        x$duration <- x$end.time - x$start.time
        # Remove any sig windows (rows) under min sig window duration
        x <- x[x$duration >= ifelse(grepl("min=50ms",alpha.loop), 50, 100),]
      }
      return(x)
    })
    elecs.with.syntax <- names(which(sapply(sig.syntax.windows, function(x){nrow(x) > 0})))
    
    
    ### Add terms ts from SP
    ## Get sample labels to keep by task
    sp.og.samples <- time.convert(rownames(zs.adjusted[[1]]), "sample.labels", "samples")
    sp.keep.these.samples <- sp.og.samples[which(sp.og.samples > (-median.rt.samples['sp']))]
    sp.keep.these.sample.labels <- time.convert(sp.keep.these.samples, "samples", "sample.labels")
    
    ## Get SP data to stack
    # diff.voice
    sp.voice <- data.frame(bind_cols(lapply(zs.adjusted, function(x){
      x[sp.keep.these.sample.labels, 'diff.voice']
    })))
    
    # Make stacked data (in parallel to stacked data for original NMF clustering)
    stacked.data[['sp']] <- sp.voice[,colnames(sp.voice)[colnames(sp.voice) %in% elecs.with.syntax]]
    
    # Clean up
    rm(sp.og.samples, sp.keep.these.samples, sp.keep.these.sample.labels, sp.voice)
    
    
    ### Combine
    stacked.data <- bind_rows(stacked.data)
    
    
    ## Select appropriate elecs depending on loops
    # All electrodes with sig elecs for this alpha loop that are also in use.these.elecs
    sp.sig.elecs.all <- intersect(unique(unlist(zs.sig.elecs[[alpha.loop]], use.names = FALSE)), use.these.elecs)
    stacked.data <- stacked.data[, colnames(stacked.data)[colnames(stacked.data) %in% sp.sig.elecs.all]]
    
    
    ##
    ## NMF
    ##
    
    message("Running NMF models! ", Sys.time())
    # Run NMF
    nmf.models <- nmf(stacked.data, 
                      rank = nmf.ranks.to.try, 
                      nrun = 50, # high - default is 30 (i.e., runs NMF 30 times per rank to find a consensus matrix)
                      .opt = 'vP', # mac (already running in parallel automatically)
                      .pbackend = n.cores.to.use.nmf,
                      seed = 123)
    
    
    # Save results
    dir.create(nmf.results.path, showWarnings = FALSE, recursive = TRUE)
    save.these <- c("nmf.models",
                    "nmf.ranks.to.try")    
    save(list = save.these,
         file = paste0(nmf.results.path, nmf.results.file))
    
  } # if(file.exists(paste0(nmf.results.path, nmf.results.file))){}else{
  
  
  ### NMF summary plots
  save.nmf.rank.summary.fig.path <- paste0(output.path,'figures/1t - clustering just syntax/',
                                           alpha.loop,'/',
                                           'all models/summary/')
  dir.create(save.nmf.rank.summary.fig.path, showWarnings = FALSE, recursive = TRUE)
  
  # Plot rank-to-rank differences
  rank.diffs <- data.frame('current.rank' = nmf.ranks.to.try[-1],
                           'previous.rank' = nmf.ranks.to.try[-length(nmf.ranks.to.try)],
                           'evar' = NA,
                           'residuals' = NA,
                           'rss' = NA,
                           'cophenetic' = NA)
  for(rank.diff.loop in 1:nrow(rank.diffs)){
    # rank.diff.loop = 1
    for(stat.loop in 3:ncol(rank.diffs)){ # stat.loop = 3
      current.stat <- colnames(rank.diffs)[stat.loop]
      rank.diffs[rank.diff.loop, stat.loop] <- 
        nmf.models$measures[nmf.models$measures$rank == rank.diffs$current.rank[rank.diff.loop], current.stat] -
        nmf.models$measures[nmf.models$measures$rank == rank.diffs$previous.rank[rank.diff.loop], current.stat]
    }; rm(stat.loop)
  }; rm(rank.diff.loop)
  pdf(paste0(save.nmf.rank.summary.fig.path,'NMF summary stats - differences from rank to rank.pdf'),
      width = 8, height = 6)
  par(mfrow = c(2,2))
  for(stat.loop in 3:ncol(rank.diffs)){
    plot(x = rank.diffs$current.rank,
         y = rank.diffs[,stat.loop],
         type = 'o',
         main = names(rank.diffs)[stat.loop],
         ylab = 'diff. from prev.',
         xlab = 'rank')
    abline(h = 0, lty = 2, col = 'grey')
  }; rm(stat.loop)
  dev.off()
  
  # Plot rank comparisons
  pdf(paste0(save.nmf.rank.summary.fig.path,'NMF summary stats.pdf'),
      width = 8,
      height = 6)
  print(plot(nmf.models))
  # # OR also try:
  # ggsave(paste0(save.nmf.rank.summary.fig.path,'NMF summary stats.pdf'), 
  #        plot = plot(nmf.models))
  dev.off()
  
  
  
  ### Cluster data for line plots and brain plots
  ## Loop thru NMFS and create data for brain plots
  for(nmf.rank.loop in nmf.ranks.to.try){
    # nmf.rank.loop = nmf.ranks.to.try[1]
    current.label = paste0('rank=',nmf.rank.loop)  
    
    # Get data from results
    factors <- data.frame(basis(nmf.models$fit[[as.character(nmf.rank.loop)]]))
    weights <- data.frame(t(coef(nmf.models$fit[[as.character(nmf.rank.loop)]])))
    colnames(factors) <- colnames(weights) <- paste0("NMF_",1:ncol(factors))
    
    # Assign elecs to clusters corresponding to the factors - whichever factor has highest value
    cluster.assignments <- data.frame('electrode' = row.names(weights),
                                      'cluster' = apply(weights, 1, which.max))
    # cluster.assignments[names(which(apply(weights, 1, max) < cluster.threshold)),]$cluster <- nmf.rank.loop + 1
    # Add unique color for each cluster
    cluster.assignments$color <- cubicl(length(unique(cluster.assignments$cluster)))[cluster.assignments$cluster]
    cluster.assignments$cluster <- paste0('NMF_',cluster.assignments$cluster)
    cluster.assignments <- cbind(cluster.assignments,
                                 t(col2rgb(cluster.assignments$color)))
    # Add unique color for each patient
    load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))
    patient.colors <- colors$white$patients[substr(cluster.assignments$electrode, 1, 5), 'hex']
    patient.colors <- data.frame(t(col2rgb(patient.colors)))
    colnames(patient.colors) <- paste0("patient_",colnames(patient.colors))
    cluster.assignments <- cbind(cluster.assignments,
                                 patient.colors)
    rm(patient.colors)
    
    # Clusters one-hot coded
    cluster.assignments.1hot <- data.frame('cluster' = cluster.assignments$cluster)
    for(cluster.loop in 1:nmf.rank.loop){ # cluster.loop = 1
      cluster.assignments.1hot[,paste0('NMF_',cluster.loop)] <-
        as.numeric(cluster.assignments.1hot$cluster == paste0('NMF_',cluster.loop))
    }; rm(cluster.loop)
    cluster.assignments.1hot$cluster <- NULL
    
    # Colors
    cluster.colors <- cubicl(nmf.rank.loop)
    names(cluster.colors) <- paste0('NMF_',1:nmf.rank.loop)
    
    # In RGB
    cluster.colors.rgb <- data.frame(cluster.colors)
    cluster.colors.rgb <- cbind(cluster.colors.rgb,
                                t(col2rgb(cluster.colors.rgb$cluster.colors)))
    cluster.colors.rgb$cluster.colors <- NULL
    cluster.colors.rgb$alpha <- 255
    cluster.colors.rgb <- cluster.colors.rgb / 255
    
    ### Make a row label matrix for localizations
    weights.and.clusters.localizations <- elec.info[row.names(weights),
                                                    c('patient',
                                                      'MNI_x','MNI_y','MNI_z',
                                                      'T1_x','T1_y','T1_z',
                                                      'region_clinical',
                                                      'active','bad_elec','bad_localization')]
    
    ### Combine
    all.cluster.plotting.data.csv <- cbind(
      cluster.assignments[,c('electrode','cluster')],
      weights.and.clusters.localizations[,c('patient','MNI_x','MNI_y','MNI_z','region_clinical')],
      weights
    )
    
    ### Add peak times and values and start times for REIs
    term.times <- all.cluster.plotting.data.csv <- cbind(
      cluster.assignments[,c('electrode','cluster')],
      weights.and.clusters.localizations[,c('patient','MNI_x','MNI_y','MNI_z','region_clinical')]
    )
    
    # For getting sig windows etc.
    current.alpha <- gsub("including_inactive_only_if_RSA_sig_","",alpha.loop)
    
    # Loop thru terms
    for(term.loop in names(zs.sig.elecs[[1]])){
      # term.loop = names(zs.sig.elecs[[1]])[1]
      
      # Get all elec-term-specific time series in dataframe
      term.vals <- list()
      for(elec.loop in term.times$elec){
        # elec.loop = term.times$elec[1]
        term.vals[[elec.loop]] <- zs.adjusted[[elec.loop]][,term.loop,drop = FALSE]
        colnames(term.vals[[elec.loop]]) <- elec.loop
      }; rm(elec.loop)
      term.vals <- bind_cols(term.vals)
      
      # Terms sig?
      term.times[,paste0(term.loop,'_sig')] <- as.numeric(term.times$electrode %in% zs.sig.elecs[[current.alpha]][[term.loop]])
      
      # Get peak RSA term value
      term.times[,paste0('peak_val_',term.loop)] <- apply(term.vals, 2, max)
      
      # Get time at RSA term value peak
      term.times[,paste0('peak_time_',term.loop)] <- 
        time.convert(rownames(term.vals)[apply(term.vals, 2, which.max)], "sample.labels", "times")
      
      # Get start time of significance for RSA term
      term.times[,paste0('start_time_',term.loop)] <- NA
      for(elec.loop in term.times$elec){
        # elec.loop = term.times$elec[1]
        
        # Get sig windows
        current.sig.windows <- 
          zs.sig.windows[[elec.loop]][[term.loop]][[current.alpha]]
        if(nrow(current.sig.windows) > 0){ # if any sig windows
          term.times[elec.loop, paste0('start_time_',term.loop)] <- min(current.sig.windows$start.time)  
        }else{ # if any sig windows
          term.times[elec.loop, paste0('start_time_',term.loop)] <- NA
        } # if any sig windows
      }; rm(elec.loop)
      
      # Get earliest -- peak or start
      term.times[,paste0('earliest_time_',term.loop)] <- NA
      for(elec.loop in term.times$elec){
        term.times[elec.loop,paste0('earliest_time_',term.loop)] <-
          min(term.times[elec.loop,paste0('start_time_',term.loop)], 
              term.times[elec.loop,paste0('peak_time_',term.loop)],
              na.rm = TRUE)
      }; rm(elec.loop)
      
    }; rm(term.loop)
    colnames(term.times) <- gsub(".","_",colnames(term.times), fixed = TRUE)
    
    ### Save results
    current.save.data.path <- paste0(output.path,'data/1t - clustering just syntax/',
                                     alpha.loop,'/',
                                     'individual models/RData files/rank=',nmf.rank.loop,'/')
    dir.create(current.save.data.path, showWarnings = FALSE, recursive = TRUE)
    save.these <- c('factors',
                    'weights',
                    'all.cluster.plotting.data.csv',
                    'cluster.assignments',
                    'cluster.assignments.1hot',
                    'cluster.colors',
                    'cluster.colors.rgb',
                    'weights.and.clusters.localizations',
                    'term.times',
                    'nmf.ranks.to.try')
    save(list = save.these,
         file = paste0(current.save.data.path,'NMF clustering results.RData'))
    
    # Save CSVs to plot brains
    current.save.data.path <- paste0(output.path,'data/1t - clustering just syntax/',
                                     alpha.loop,'/',
                                     'individual models/csvs for brain plots/rank=',nmf.rank.loop,'/')
    dir.create(current.save.data.path, showWarnings = FALSE, recursive = TRUE)
    write.csv(weights,
              file = paste0(current.save.data.path, 'weights.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(cluster.assignments,
              file = paste0(current.save.data.path, 'cluster_assignments.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(cluster.assignments.1hot,
              file = paste0(current.save.data.path, 'cluster_assignments_1hot.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(weights.and.clusters.localizations,
              file = paste0(current.save.data.path, 'row_labels_and_localizations.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(cluster.colors.rgb,
              file = paste0(current.save.data.path, 'cluster_colors.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(all.cluster.plotting.data.csv,
              file = paste0(current.save.data.path, 'all_cluster_plotting_data.csv'),
              quote = FALSE, row.names = FALSE)
    write.csv(term.times,
              file = paste0(current.save.data.path, 'term_times.csv'),
              quote = FALSE, row.names = FALSE)
    
    
    ### Quick cluster plots
    term.colors <- c(setNames(colors$black$rsa_term_line_plots$hex, 
                              rownames(colors$black$rsa_term_line_plots)),
                     c('diff.log.rt' = colors$white$rainbow_bright['grey','hex']))
    plot.grid <- c(2, nmf.rank.loop)
    save.grid.plot.path <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                  alpha.loop,'/',
                                  'individual models/quick and dirty grid plots/')
    dir.create(save.grid.plot.path, showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(save.grid.plot.path, 'rank=', nmf.rank.loop,'.pdf'),
        height = 3 * plot.grid[1],
        width = 4 * plot.grid[2])
    par(mfcol = plot.grid,
        oma = c(0,0,1,0),
        bg = 'black')
    for(cluster.loop in sort(unique(cluster.assignments$cluster))){
      # cluster.loop = unique(cluster.assignments$cluster)[1]
      
      current.elecs <- lapply(split(cluster.assignments, cluster.assignments$cluster), rownames)[[cluster.loop]]
      
      # ECoG in current band.loop
      current.elec.means <- lapply(elec.means, function(x){apply(x[,current.elecs], 1, mean, na.rm = TRUE)})
      current.elec.ses <- lapply(elec.means, function(x){apply(x[,current.elecs], 1, function(y){sd(y, na.rm = TRUE) / sqrt(length(y[!is.na(y)]))})})
      current.xs <- lapply(current.elec.means, names)
      plot.time.series(.y.values = current.elec.means,
                       .x.values = current.xs,
                       .error.bars = current.elec.ses,
                       .x.limits = c(time.convert(-max(median.rt.samples), "samples", "times") - 100, 300),
                       .y.limits.min.at.least = ifelse(band.loop == 'high_gamma', -1, -.25),
                       .y.limits.max.at.least = ifelse(band.loop == 'high_gamma', 3, 1.25),
                       .y.ticks = 0:ifelse(band.loop == 'high_gamma',3,1),
                       .colors = colors$white$task_line_plots[names(current.elec.means),'hex'],
                       .horizontal.line.at = 0,
                       .y.label = ifelse(band.loop == 'beta', 'Beta','HGA'),
                       .title = paste0(cluster.loop, " (", length(current.elecs), " elecs)"))
      
      # Z-scored t-values
      current.z.means <- as.list(data.frame(elementwise.matrix.apply(zs.adjusted[current.elecs], .function = "mean")))
      current.z.xs <- rownames(zs.adjusted[[1]])
      plot.time.series(.y.values = current.z.means,
                       .x.values = current.z.xs,
                       .x.limits = c(time.convert(-max(median.rt.samples), "samples", "times") - 100, 300),
                       .y.limits.min.at.least = -1,
                       .y.limits.max.at.least = 3,
                       .colors = term.colors[names(current.z.means)],
                       .horizontal.line.at = 0,
                       .y.label = 'importance')
      
      # # ECoG in other band.loop
      # current.elec.means <- lapply(other.band.elec.means, function(x){apply(x[,current.elecs], 1, mean, na.rm = TRUE)})
      # current.elec.ses <- lapply(other.band.elec.means, function(x){apply(x[,current.elecs], 1, function(y){sd(y, na.rm = TRUE) / sqrt(length(y[!is.na(y)]))})})
      # current.xs <- lapply(current.elec.means, names)
      # plot.time.series(.y.values = current.elec.means,
      #                  .x.values = current.xs,
      #                  .error.bars = current.elec.ses,
      #                  .x.limits = c(time.convert(-max(median.rt.samples), "samples", "times") - 100, 300),
      #                  .y.limits.min.at.least = ifelse(band.loop == 'high_gamma', -.25, -1),
      #                  .y.limits.max.at.least = ifelse(band.loop == 'high_gamma', 1, 3),
      #                  .y.ticks = 0:ifelse(band.loop == 'high_gamma',1,3),
      #                  .colors = colors$white$task_line_plots[names(current.elec.means),'hex'],
      #                  .horizontal.line.at = 0,
      #                  .y.label = ifelse(band.loop == 'beta', 'HGA','Beta'))
      
    }; rm(cluster.loop)
    dev.off()
    
    
    ### Plot individual clusters - high gamma AND RSA terms together
    term.colors <- c(setNames(colors$white$rsa_term_line_plots$hex, 
                              rownames(colors$white$rsa_term_line_plots)),
                     c('diff.log.rt' = colors$white$rainbow_bright['grey','hex']))
    individual.plot.grid <- c(2, 1)
    save.individual.plot.path <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                        alpha.loop,'/',
                                        'individual models/quick and dirty individual plots/',
                                        'rank=',nmf.rank.loop,'/')
    dir.create(save.individual.plot.path, showWarnings = FALSE, recursive = TRUE)
    for(cluster.loop in sort(unique(cluster.assignments$cluster))){
      # cluster.loop = unique(cluster.assignments$cluster)[1]
      pdf(paste0(save.individual.plot.path, 'rank=', nmf.rank.loop,' - cluster=',cluster.loop,'.pdf'),
          height = 3 * individual.plot.grid[1],
          width = 4 * individual.plot.grid[2])
      par(mfcol = individual.plot.grid,
          oma = c(0,0,1,0),
          bg = 'white')
      
      current.elecs <- lapply(split(cluster.assignments, cluster.assignments$cluster), rownames)[[cluster.loop]]
      
      # ECoG in current band.loop
      current.elec.means <- lapply(elec.means, function(x){apply(x[,current.elecs], 1, mean, na.rm = TRUE)})
      current.elec.ses <- lapply(elec.means, function(x){apply(x[,current.elecs], 1, function(y){sd(y, na.rm = TRUE) / sqrt(length(y[!is.na(y)]))})})
      current.xs <- lapply(current.elec.means, names)
      plot.time.series(.y.values = current.elec.means,
                       .x.values = current.xs,
                       .error.bars = current.elec.ses,
                       .x.limits = c(time.convert(-max(median.rt.samples), "samples", "times") - 100, 300),
                       .y.limits.min.at.least = ifelse(band.loop == 'high_gamma', -1, -.25),
                       .y.limits.max.at.least = ifelse(band.loop == 'high_gamma', 3, 1.25),
                       .y.ticks = 0:ifelse(band.loop == 'high_gamma',3,1),
                       .colors = colors$white$task_line_plots[names(current.elec.means),'hex'],
                       .theme = 'white',
                       .horizontal.line.at = 0,
                       .y.label = ifelse(band.loop == 'beta', 'Beta','HGA'),
                       .title = paste0(cluster.loop, " (", length(current.elecs), " elecs)"))
      
      # Z-scored t-values
      current.z.means <- as.list(data.frame(elementwise.matrix.apply(zs.adjusted[current.elecs], .function = "mean")))
      current.z.xs <- rownames(zs.adjusted[[1]])
      plot.time.series(.y.values = current.z.means,
                       .x.values = current.z.xs,
                       .x.limits = c(time.convert(-max(median.rt.samples), "samples", "times") - 100, 300),
                       .y.limits.min.at.least = -1,
                       .y.limits.max.at.least = 3,
                       .colors = term.colors[names(current.z.means)],
                       .theme = 'white',
                       .horizontal.line.at = 0,
                       .y.label = 'importance')
      
      dev.off()
    }; rm(cluster.loop)
    
    
  }; rm(nmf.rank.loop)
  
  # Clean up
  rm(list = save.these)
  gc()
  
} # alpha.loop
# }; rm(band.loop)



# Finish!
message('Script completed successfully. ',Sys.time())











