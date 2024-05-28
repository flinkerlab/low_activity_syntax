### Correlate RSIs and high gamma across electrodes at each time sample
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
    n.cores.to.use = 18
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
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Get min and max
source(paste0(path,'/analysis/R/functions/min_max.R'))
# Violin plot
source(paste0(path,'/analysis/R/functions/plot_violin.R'))


### Colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

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
               (elec.info$bad_elec == 0) &
               (elec.info$visual_elec == 0) &
               (elec.info$active == 1))
rownames(elec.info) <- elec.info$patient_elec
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec

## Add region supercategories
gross.rois.df <- read.csv('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/define ROIs/output/gross ROIs.csv')
rownames(gross.rois.df) <- gross.rois.df$region.clinical
elec.info$gross.roi <- gross.rois.df[elec.info$region_clinical,'region.category']


### Loop thru beta/high gamma data
# for(band.loop in c('high_gamma','beta')[1]){ ## UNCOMMENT
band.loop = c('high_gamma','beta')[1] ## UNCOMMENT


# Directories
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/')


##
## Load all data
##

### Get median.rt.samples
load(paste0(path, 'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))

### Elec means - high gamma
load(paste0(path,'/analysis/R/electrode means/means - warped data - multiband/output - high_gamma/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses

# Get metadata
tasks <- names(elec.means)

# Add NA dataframe for Patient006 LP block
blank.Patient006.lp.df <- data.frame(matrix(nrow = nrow(elec.means$lp[[1]]),
                                       ncol = ncol(elec.means$pn[['Patient006']])))
colnames(blank.Patient006.lp.df) <- colnames(elec.means$pn[['Patient006']])
rownames(blank.Patient006.lp.df) <- rownames(elec.means$lp[[1]])
elec.means$lp[['Patient006']] <- elec.sds$lp[['Patient006']] <- elec.ses$lp[['Patient006']] <- blank.Patient006.lp.df
rm(blank.Patient006.lp.df)

# Collapse across patients
elec.means <- lapply(elec.means, bind_cols)
# elec.sds <- lapply(elec.sds, bind_cols)
# elec.ses <- lapply(elec.ses, bind_cols)

# Smooth to match smoothing profile of z-values (40 samples x 2)
elec.means <- lapply(elec.means, function(x){x <- data.frame(apply(x, 2, smoothing, n.samples.pre = 40))})
# elec.sds <- lapply(elec.sds, function(x){x <- data.frame(apply(x, 2, smoothing, n.samples.pre = 40))})
# elec.ses <- lapply(elec.ses, function(x){x <- data.frame(apply(x, 2, smoothing, n.samples.pre = 40))})

### Load unadjusted smoothed z-values
unadjusted.ts.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
unadjusted.ts.file <- 'just z-scores.RData'
load(paste0(unadjusted.ts.path, unadjusted.ts.file)) # loads ts.z

# Terms
terms <- names(ts.z[[1]])
term.labels <- c('diff.voice' = 'syntax',
                 'diff.1st.word' = 'word',
                 'diff.event.semantics' = 'semantics')

# Set negative unadjusted z-values to 0
nonneg.ts.z <- lapply(ts.z, function(x){
  data.frame(apply(x, 2, function(y){
    y[y<0] <- 0
    return(y)}))})

# If any elec term is all 0s, replace with NAs so it doesn't affect correlations
nonneg.ts.z <- 
  lapply(nonneg.ts.z, function(x){
    for(term.loop in colnames(x)){
      if(all(x[,term.loop] == 0)){
        x[,term.loop] <- NA
      } # if(all(x[,term.loop] == 0)){
    }; rm(term.loop)
    return(x)
  })

### Load adjusted smoothed z-values
adjusted.ts.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
adjusted.ts.file <- 'adjusted zs with sig windows and elecs.RData'
load(paste0(adjusted.ts.path, adjusted.ts.file))

# Sample labels to plot for each task
keep.sample.labels <- list()
for(term.loop in names(median.rt.samples)){
  last.sample <- time.convert(500, "times", "samples")
  keep.sample.labels[[term.loop]] <- time.convert((-median.rt.samples[term.loop] - 100):last.sample,
                                                  "samples", "sample.labels")
  rm(last.sample)
}; rm(term.loop)

# Common sample labels
common.sample.labels <- intersect(rownames(nonneg.ts.z[[1]]), rownames(elec.means$pn))

# Elecs present in all three datasets
elecs.to.loop.thru <- intersect(names(nonneg.ts.z), names(elec.means$pn))

# Metadata
terms <- colnames(zs.adjusted[[1]])
term.labels <- c('diff.voice' = 'syntax',
                 'diff.event.semantics' = 'semantics',
                 'diff.1st.word' = 'word')

### Loop through significance thresholds for elecs
sig.elec.hps <- c('alpha=.05_min=100ms',
                  'alpha=.01_min=50ms')

keep <- c(ls(), "keep", "sig.hp.loop")

  
  # for(sig.hp.loop in sig.elec.hps){ ## UNCOMMENT
  sig.hp.loop = sig.elec.hps[1] ## UNCOMMENT
  
  # Clean up
  rm(list = ls()[which(!ls() %in% keep)])
  gc()
  
  ### Get correlation of ECoG activity and RSA terms for each elec
  # Storage
  ecog.rsa.cors <- 
    ecog.rsa.sig <- 
    data.frame(matrix(nrow = length(common.sample.labels),
                      ncol = length(terms),
                      dimnames = list(common.sample.labels, terms)))
  ecog.rsa.sig.windows <- list()
  
  # Loop!
  for(term.loop in terms){
    # term.loop = terms[1]
    
    # Get sig elecs this term
    current.elecs <- zs.sig.elecs[[sig.hp.loop]][[term.loop]]
    
    for(sample.loop in common.sample.labels){
      # sample.loop = common.sample.labels[1]
      
      # Get RSA and ECoG vals at this sample for the sig elecs for this term
      current.zs <- sapply(zs.adjusted[current.elecs], function(x){
        return(x[sample.loop, term.loop])})
      current.ecog <- unlist(elec.means$sp[sample.loop, current.elecs])
      
      # Correlate
      current.cor <- cor.test(current.zs, current.ecog, 
                              method = 'spearman')
      ecog.rsa.cors[sample.loop, term.loop] <- current.cor$estimate
      ecog.rsa.sig[sample.loop, term.loop] <- 
        as.numeric((current.cor$p.value < .05) & # sig and
                     (current.cor$estimate > 0)) # positive relationship
      
    }#; rm(sample.loop)  
    
    ecog.rsa.sig.windows[[term.loop]] <- 
      get.significant.windows(ecog.rsa.sig[,term.loop],
                              .sample.labels = rownames(ecog.rsa.sig),
                              output.class = "data.frame",
                              include.duration = TRUE,
                              .exclude.sig.durations.under.ms = 100,
                              .exclude.times.before.ms = time.convert(-median.rt.samples['sp'], "samples", "times"))
  }#; rm(term.loop)  
  
  
  
  ###
  ### Plot results
  ###
  
  current.x.min <- min(time.convert(common.sample.labels, "sample.labels", "times"))
  current.x.max <- current.x.min + length(-1600:1100) # width of plots from "/task ECoG comparisons/ROI wilcox tests and squiggle plots/1a - squiggle plots - regions - data unwarped.R"
  
  
  # for(theme.loop in c('white','black')){ ## UNCOMMENT
  theme.loop = c('white','black')[1] ## UNCOMMENT
  
  # for(transparency.loop in c(TRUE, FALSE)){ ## UNCOMMENT
  transparency.loop = c(TRUE, FALSE)[1] ## UNCOMMENT
  
  ##
  ## PDF for visualization - labels etc.
  ## 
  
  save.dir <- paste0(output.path, 
                     'figures/1k - ECoG RSA relationship/1k_c - ECoG-RSA cross-elec correlations/',
                     'with hga/',
                     'time series/visualization/',
                     theme.loop,' - ',ifelse(transparency.loop, 'transparent/', 'opaque/'))
  dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
  
  ### Correlation r
  pdf(paste0(save.dir, 'ECoG-RSA cross-elec correlation r - ',sig.hp.loop,'.pdf'),
      width = 7, height = 5)
  par(oma = c(1,0,1,0) * 1.8)
  plot.time.series(.y.values = as.list(ecog.rsa.cors),
                   .x.values = rownames(ecog.rsa.cors),
                   show.t0 = TRUE,
                   .horizontal.line.at = 0,
                   .colors = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .sig.windows = ecog.rsa.sig.windows[names(as.list(ecog.rsa.cors))],
                   .sig.color = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .y.label = 'correlation (r)',
                   .x.label = 'time (ms)',
                   .y.limits.min.at.least = -.6,
                   .y.limits.max.at.least = .6,
                   .y.ticks = c(-.5, 0, .5),
                   .x.limits = c(current.x.min, 500),
                   .x.ticks = c(-1000, -500, 0, 500),
                   .x.tick.labels = c('-1000', '', '0', '500'),
                   .zoom = 1.8,
                   .theme = theme.loop,
                   .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
  add.text.line.multiple.colors(term.labels[names(as.list(ecog.rsa.cors))],
                                colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                                .outer = TRUE)
  dev.off()
  
  ### R^2
  current.ys <- as.list(ecog.rsa.cors^2)
  current.max <- max(unlist(current.ys))
  current.order.of.mag <- 10
  current.y.upper.lim <- NA
  while(is.na(current.y.upper.lim)){
    if(floor(current.max * current.order.of.mag) > 0){
      current.y.upper.lim <- ceiling(current.max * current.order.of.mag) / current.order.of.mag
    }else{ # if(floor(current.max * current.order.of.mag) > 0){
      current.order.of.mag <- current.order.of.mag * 10
    } # if(floor(current.max * current.order.of.mag) > 0){}else{
  } # while(is.na(current.y.upper.lim)){
  pdf(paste0(save.dir, 'ECoG-RSA cross-elec correlation r-sqaured - ',sig.hp.loop,'.pdf'),
      width = 7, height = 5)
  par(oma = c(1,0,1,0) * 1.8)
  plot.time.series(.y.values = as.list(ecog.rsa.cors^2),
                   .x.values = rownames(ecog.rsa.cors),
                   show.t0 = TRUE,
                   .horizontal.line.at = 0,
                   .colors = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .sig.windows = ecog.rsa.sig.windows[names(as.list(ecog.rsa.cors))],
                   .sig.color = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .y.label = 'correlation (r^2)',
                   .x.label = 'time (ms)',
                   .x.limits = c(current.x.min, 500),
                   .y.limits = c(0, current.y.upper.lim),
                   .y.ticks = c(0, current.y.upper.lim),
                   .x.ticks = c(-1000, -500, 0, 500),
                   .x.tick.labels = c('-1000', '', '0', '500'),
                   .zoom = 1.8,
                   .theme = theme.loop,
                   .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
  add.text.line.multiple.colors(term.labels[names(as.list(ecog.rsa.cors))],
                                colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                                .outer = TRUE)
  dev.off()
  
  
  ##
  ## PDF for pub - no axis labels etc.
  ## 
  
  save.dir <- paste0(output.path, 
                     'figures/1k - ECoG RSA relationship/1k_c - ECoG-RSA cross-elec correlations/',
                     'with hga/',
                     'time series/PDF for publication/',
                     theme.loop,' - ',ifelse(transparency.loop, 'transparent/', 'opaque/'))
  save.dir.terms <- paste0(output.path, 
                     'figures/1k - ECoG RSA relationship/1k_c - ECoG-RSA cross-elec correlations/',
                     'with hga/',
                     'time series/PDF for publication - individual terms/',
                     theme.loop,' - ',ifelse(transparency.loop, 'transparent/', 'opaque/'))
  dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(save.dir.terms, showWarnings = FALSE, recursive = TRUE)
  
  ### Correlation r
  pdf(paste0(save.dir, 'ECoG-RSA cross-elec correlation r - ',sig.hp.loop,'.pdf'),
      width = 6, height = 4)
  par(oma = c(0,0,1,0))
  plot.time.series(.y.values = as.list(ecog.rsa.cors),
                   .x.values = rownames(ecog.rsa.cors),
                   show.t0 = TRUE,
                   .horizontal.line.at = 0,
                   .colors = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .sig.windows = ecog.rsa.sig.windows[names(as.list(ecog.rsa.cors))],
                   .sig.color = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .y.label = '',
                   .x.label = '',
                   .y.limits = c(-1, 1),
                   .y.ticks = c(-1, 0, 1),
                   # .x.limits = c(current.x.min, current.x.max),
                   .x.limits = c(time.convert(-median.rt.samples['sp'], "samples", "times") - 150, 500),
                   .x.ticks = c(-1000, -500, 0, 500),
                   .x.tick.labels = c('-1000', '', '0', '500'),
                   .zoom = 1.8,
                   .margin = c(3,3,0,1),
                   .theme = theme.loop,
                   .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
  dev.off()
  
  ### Correlation r
  for(width.loop in c(6, 10)){
    # width.loop = 6
    pdf(paste0(save.dir, 'ECoG-RSA cross-elec correlation r - custom y-limits - width = ',width.loop,' - ',sig.hp.loop,'.pdf'),
        width = width.loop, height = 4.5) # 5.5 (was too tall try 5.3)
    par(oma = c(0,0,2,0))
    y.limits <- max(.5, max(unlist(ecog.rsa.cors))) * c(-1, 1)
    y.ticks <- max(floor(10*y.limits[2]),1)/10 * c(-1,0,1)
    plot.time.series(.y.values = as.list(ecog.rsa.cors),
                     .x.values = rownames(ecog.rsa.cors),
                     show.t0 = TRUE,
                     .horizontal.line.at = 0,
                     .colors = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                     .sig.windows = ecog.rsa.sig.windows[names(as.list(ecog.rsa.cors))],
                     .sig.color = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                     .y.label = '',
                     .x.label = '',
                     .y.limits = y.limits,
                     .y.ticks = y.ticks,
                     # .x.limits = c(current.x.min, current.x.max),
                     .x.limits = c(time.convert(-median.rt.samples['sp'], "samples", "times") - 150, 500),
                     .x.ticks = c(-1000, -500, 0, 500),
                     # .x.tick.labels = c('-1000', '', '0', '500'),
                     .zoom = 1.8,
                     .margin = c(3,3,0,1),
                     .theme = theme.loop,
                     .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
    dev.off()
    
    # By term
    for(term.loop in colnames(ecog.rsa.cors)){
      pdf(paste0(save.dir.terms, 'ECoG-RSA cross-elec corr r_custom y_width=',width.loop,'_',sig.hp.loop,'_',term.loop,'.pdf'),
          width = width.loop, height = 4.3) # 5.5 (was too tall try 5.3)
      par(oma = c(0,0,2,0))
      y.limits <- max(.5, max(unlist(ecog.rsa.cors))) * c(-1, 1)
      y.ticks <- max(floor(10*y.limits[2]),1)/10 * c(-1,0,1)
      plot.time.series(.y.values = as.list(ecog.rsa.cors)[[term.loop]],
                       .x.values = rownames(ecog.rsa.cors),
                       show.t0 = TRUE,
                       .horizontal.line.at = 0,
                       .colors = colors[[theme.loop]]$rsa_term_line_plots[term.loop,'hex'],
                       .sig.windows = ecog.rsa.sig.windows[term.loop],
                       .sig.color = colors[[theme.loop]]$rsa_term_line_plots[term.loop,'hex'],
                       .y.label = '',
                       .x.label = '',
                       .y.limits = y.limits,
                       .y.ticks = y.ticks,
                       # .x.limits = c(current.x.min, current.x.max),
                       .x.limits = c(time.convert(-median.rt.samples['sp'], "samples", "times") - 150, 500),
                       .x.ticks = c(-1000, -500, 0, 500),
                       # .x.tick.labels = c('-1000', '', '0', '500'),
                       .zoom = 1.8,
                       .margin = c(3,3,0,1),
                       .theme = theme.loop,
                       .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
      dev.off()
    }
  }#; rm(width.loop)
  
  ### Correlation r^2
  pdf(paste0(save.dir, 'ECoG-RSA cross-elec correlation r-squared - ',sig.hp.loop,'.pdf'),
      width = 8, height = 4)
  par(oma = c(0,0,0,0))
  plot.time.series(.y.values = as.list(ecog.rsa.cors^2),
                   .x.values = rownames(ecog.rsa.cors),
                   show.t0 = TRUE,
                   .horizontal.line.at = 0,
                   .colors = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .sig.windows = ecog.rsa.sig.windows[names(as.list(ecog.rsa.cors))],
                   .sig.color = colors[[theme.loop]]$rsa_term_line_plots[names(as.list(ecog.rsa.cors)),'hex'],
                   .y.label = '',
                   .x.label = '',
                   .y.limits = c(0, current.y.upper.lim),
                   .y.ticks = c(0, current.y.upper.lim),
                   # .x.limits = c(current.x.min, current.x.max),
                   .x.limits = c(time.convert(-median.rt.samples['sp'], "samples", "times") - 150, 500),
                   .x.ticks = c(-1000, -500, 0, 500),
                   .x.tick.labels = c('-1000', '', '0', '500'),
                   .zoom = 1.8,
                   .margin = c(3,3,.5,.5),
                   .theme = theme.loop,
                   .background = ifelse(transparency.loop, rgb(1,1,1,0), theme.loop))
  dev.off()
  
  
}#; rm(transparency.loop)

}#; rm(theme.loop)

}#; rm(sig.hp.loop)

# }; rm(band.loop)

# Finish!
message('Script completed successfully. ',Sys.time())







