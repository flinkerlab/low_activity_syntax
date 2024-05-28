### Compare tasks (e.g., sentence and list production) by ROI using unwarped data
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
library('viridis') # for viridis color schemes
library('pals') # for ocean.balance, ocean.delta, ocean.curl colors
#pal.bands(ocean.balance, ocean.delta, ocean.curl)
library('scico') # for palettes cork, vik, lisbon, tofino, berlin, roma
library('reticulate') # for Python
library('stringr') # for str_split()
library('lme4') # for lmer()
library('lmerTest') # for lmer() that gives p-value
library('NMF') # for non-negative matrix factorization

# Clean up
rm(list=ls())
cat("\014")
message("Begin plotting ECoG for each patient, electrode, task, and time-lock. ",Sys.time())

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
    n.cores.to.use = 18
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
# Adjust color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add a line of text to plot with varying font color
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Get significant windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Good trial subset functions
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))


### Set up elec info
# Read in elecs
elec.info <- 
  read.csv(paste0(path,
                  '/analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
# Subset to just good elecs:
elec.info <- droplevels(subset(elec.info,
                               (bad_elec == 0) &
                                 (active.stim.locked == 1) &
                                 (active.prod.locked == 1) &
                                 (bad_localization == 0) &
                                 (! region_clinical %in% c('NaN','Unknown',''))))
use.these.elecs <- elec.info$patient_elec

# Define regions to look at
attach(paste0(path,'analysis/R/define ROIs/output/all region definitions.RData'))
region.categories <- main.rois[c('IFG','MFG','STG','MTG','IPL','pericentral')]
detach()

# Get rows in data for each of these regions
region.rows <- list()
region.elecs <- list()
for(i in 1:length(region.categories)){
  current.region <- names(region.categories)[i]
  region.rows[[current.region]] <- which(elec.info$region_clinical %in% region.categories[[current.region]])
  region.elecs[[current.region]] <- elec.info$patient_elec[region.rows[[current.region]]]
  rm(current.region)
}; rm(i)
# Turn elecs into a dataframe for easy lookup
region.elecs.df <- data.frame('elec' = unlist(region.elecs))
region.elecs.df$region.category <- rownames(region.elecs.df)
for(i in 0:9){region.elecs.df$region.category <- gsub(i, '', region.elecs.df$region.category)}
row.names(region.elecs.df) <- NULL

## Get rid of regions with fewer than 3 elecs
n.elecs.per.region <- sapply(region.elecs, length)
remove.regions <- names(n.elecs.per.region[which(n.elecs.per.region < 3)])
region.elecs.df <- droplevels(subset(region.elecs.df, !region.category %in% remove.regions))
region.elecs[remove.regions] <- NULL
region.rows[remove.regions] <- NULL
region.categories[remove.regions] <- NULL

### Metadata
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
tasks <- c('pn','sp','lp') # picture naming, sentence production, list production

# Colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

keep <- c(ls(), 'keep', 'band.loop')
for(band.loop in c('high_gamma','beta')){
  # band.loop = c('high_gamma','beta')[2]
  
  rm(list = ls()[which(! ls() %in% keep)])
  
  ### Begin reading in stats and stuff
  ## Takes ~2 hours so if possible load previous results
  save.data.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/data/',band.loop,'/region stats/')
  if(file.exists(paste0(save.data.dir, 'ecog wilcoxon task tests.RData'))){
    message('Not running stats; loading from previous save.')

    # Load data
    load(paste0(save.data.dir, 'ecog wilcoxon task tests.RData'))

  }else{
    
    # Initialize storage of stats etc
    region.means <- list()
    region.ses <- list()
    sentence.vs.list.sig.windows <- list()
    sentence.vs.naming.sig.windows <- list()
    list.vs.naming.sig.windows <- list()
    
      message("Begin data set up for loops: ",band.loop,'. ', Sys.time())
      
      ## Read in electrode data
      elec.data <- list()
      trial.info <- list()
      message('Attaching warped data...')
      attach(paste0(path,
                  'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',
                  band.loop,
                  '/data/elec data with RTs stretched to global median by task - bad trials and .025 to .95 RT outliers excluded.RData')
             ) # Loads: median.rt.samples, patient.trial.info, stretch.data, stretch.sample.labels
      for(task.loop in tasks){
        elec.data[[task.loop]] <- unlist(lapply(stretch.data, function(x){x[[task.loop]]}), recursive = FALSE)
        names(elec.data[[task.loop]]) <- gsub("^[^.]*\\.", "", names(elec.data[[task.loop]]))
        trial.info[[task.loop]] <- lapply(patient.trial.info, function(x){x[[task.loop]]})
      }; rm(task.loop)
      median.rt.samples <- median.rt.samples
      stretch.sample.labels <- stretch.sample.labels
      detach()
      message('...detached!')
      
      # Subset to just time samples you care about
      sample.labels <- list('locked_to_production_onset' = list(),
                            'locked_to_stimulus_onset' = list())
      for(task.loop in tasks){
        # task.loop = tasks[1]
        
        prod.locked.start.sample <- (-median.rt.samples[task.loop] - time.convert(100, "times", "samples"))
        prod.locked.end.sample <- time.convert(1100, "times", "samples")
        sample.labels$locked_to_production_onset[[task.loop]] <- 
          time.convert(prod.locked.start.sample:prod.locked.end.sample, "samples", "sample.labels")
        sample.labels$locked_to_stimulus_onset[[task.loop]] <- 
          stretch.sample.labels[[task.loop]]$locked_to_stimulus_onset[
            which(stretch.sample.labels[[task.loop]]$locked_to_production_onset %in% 
                sample.labels$locked_to_production_onset[[task.loop]])]
        
        elec.data[[task.loop]] <- lapply(elec.data[[task.loop]], function(x){x[, sample.labels[['locked_to_production_onset']][[task.loop]]]})
      }; rm(task.loop)
      
      # Get metadata
      times <- lapply(sample.labels, function(x){lapply(x, function(y){time.convert(y, 'sample.labels', 'times')})})
      
      
      ###
      ### Get ECoG means and SEs by region
      ###
      
      # Initialize storage for ECoG means
      region.elec.data <- list()
      region.elec.means <- list()
      region.means <- list()
      region.ses <- list()
      
      # Loop thru tasks
      for(task.loop in tasks){
        # task.loop = tasks[3]
        region.elec.data[[task.loop]] <- list()
        region.elec.means[[task.loop]] <- list()
        region.means[[task.loop]] <- list()
        region.ses[[task.loop]] <- list()
        
        # Loop thru region categories
        for(region.loop in 1:length(region.categories)){
          # region.loop = 1
          current.region <- names(region.categories)[region.loop]
          current.elecs <- region.elecs.df$elec[region.elecs.df$region.category == current.region]
          # Get rid of elecs with no data for this task
          current.elecs <- current.elecs[current.elecs %in% names(elec.data[[task.loop]])]
          
          # Store data for stats below
          region.elec.data[[task.loop]][[current.region]] <- elec.data[[task.loop]][current.elecs]
          current.n.rows <- sapply(region.elec.data[[task.loop]][[current.region]], nrow)
          # Add electrode and task to each elec dataframe
          for(elec.loop in current.elecs){
            # elec.loop <- current.elecs[1]
            region.elec.data[[task.loop]][[current.region]][[elec.loop]] <- 
              cbind(data.frame('elec' = rep(elec.loop, times = current.n.rows[elec.loop]),
                               'task' = rep(task.loop, times = current.n.rows[elec.loop])),
                    region.elec.data[[task.loop]][[current.region]][[elec.loop]])
          }; rm(elec.loop)
          # Collapse elec dataframes for this region into one big one
          region.elec.data[[task.loop]][[current.region]] <-
            do.call(rbind, region.elec.data[[task.loop]][[current.region]])
          
          # Get means for each electrode
          region.elec.means[[task.loop]][[current.region]] <- 
            data.frame(bind_rows(lapply(elec.data[[task.loop]][current.elecs],
                                        function(x){apply(x, 2, mean)})))
          # Get rid of electrodes where there's no data
          region.elec.means[[task.loop]][[current.region]] <- 
            na.omit(region.elec.means[[task.loop]][[current.region]])
          # SE across elecs
          region.ses[[task.loop]][[current.region]] <- 
            apply(region.elec.means[[task.loop]][[current.region]], 2, sd) /
            sqrt(nrow(region.elec.means[[task.loop]][[current.region]]))
          # Average of elecs
          region.means[[task.loop]][[current.region]] <- 
            apply(region.elec.means[[task.loop]][[current.region]], 2, mean)
        }; rm(region.loop)
      }; rm(task.loop)
      
      
      ### Get stats on differences between tasks
      for(lock.loop in names(sample.labels)){
        # lock.loop = names(sample.labels)[1]
      
      # Storage
      sentence.vs.list <- list()
      sentence.vs.naming <- list()
      list.vs.naming <- list()
      sentence.vs.list.sig.windows[[lock.loop]] <- list()
      sentence.vs.naming.sig.windows[[lock.loop]] <- list()
      list.vs.naming.sig.windows[[lock.loop]] <- list()
      
      # Loop thru regions
      for(region.loop in names(region.categories)){
        # region.loop = 'MTG'#names(region.categories)[1]
        message("Beginning EGoG activity models for ",region.loop," (",which(names(region.categories) == region.loop)," of ",length(region.categories)," regions), ",lock.loop,". ",Sys.time())
        
        # Storage
        sentence.vs.list[[region.loop]] <- list()
        sentence.vs.naming[[region.loop]] <- list()
        list.vs.naming[[region.loop]] <- list()
        sentence.vs.list.sig.windows[[lock.loop]][[region.loop]] <- list()
        sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
        list.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
        
        # Loop thru time samples
        for(sample.loop in sample.labels[[lock.loop]][['sp']]){ # sp has the most sample labels
          # sample.loop = sample.labels[[lock.loop]][['sp']][1]
          
          if(lock.loop == 'locked_to_production_onset'){
            # If locked to production onset
            current.prod.locked.sample.sp <- sample.loop
            current.prod.locked.sample.lp <- sample.labels[['locked_to_production_onset']][['lp']][
              which(sample.labels[['locked_to_production_onset']][['lp']] == sample.loop)]
            current.prod.locked.sample.pn <- sample.labels[['locked_to_production_onset']][['pn']][
              which(sample.labels[['locked_to_production_onset']][['pn']] == sample.loop)]
          }else{ # if(lock.loop == 'locked_to_production_onset'){
            # If locked to stimulus onset
            current.prod.locked.sample.sp <- 
              sample.labels[['locked_to_production_onset']][['sp']][
                which(sample.labels[['locked_to_stimulus_onset']][['sp']] == sample.loop)]
            current.prod.locked.sample.lp <- 
              sample.labels[['locked_to_production_onset']][['lp']][
                which(sample.labels[['locked_to_stimulus_onset']][['lp']] == sample.loop)]
            current.prod.locked.sample.pn <- 
              sample.labels[['locked_to_production_onset']][['pn']][
                which(sample.labels[['locked_to_stimulus_onset']][['pn']] == sample.loop)]
          } # if(lock.loop == 'locked_to_production_onset'){
          
          ### Sentence vs. list
          if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.lp)) == 2){
            current.model <-
              wilcox.test(region.elec.data$sp[[region.loop]][, current.prod.locked.sample.sp],
                          region.elec.data$lp[[region.loop]][, current.prod.locked.sample.lp],
                          alternative = 'greater')
            
            # Store p-value
            sentence.vs.list[[region.loop]][[sample.loop]] <-
              current.model$p.value
            rm(current.model)
          } # if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.lp)) == 2){
          
          
          ### Sentence vs. naming
          if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.pn)) == 2){
            current.model <- 
              wilcox.test(region.elec.data$sp[[region.loop]][, current.prod.locked.sample.sp],
                          region.elec.data$pn[[region.loop]][, current.prod.locked.sample.pn],
                          alternative = 'greater')
            
            # Store p-value
            sentence.vs.naming[[region.loop]][[sample.loop]] <- 
              current.model$p.value
            rm(current.model)
          } # if(sample.loop %in% sample.labels[[lock.loop]][['pn']]){
          
          ### List vs. naming
          if(length(c(current.prod.locked.sample.lp, current.prod.locked.sample.pn)) == 2){
            current.model <- 
              wilcox.test(region.elec.data$lp[[region.loop]][, current.prod.locked.sample.lp],
                          region.elec.data$pn[[region.loop]][, current.prod.locked.sample.pn],
                          alternative = 'greater')
            
            # Store p-value
            list.vs.naming[[region.loop]][[sample.loop]] <-
              current.model$p.value
            rm(current.model)
          } # if(sample.loop %in% sample.labels[[lock.loop]][['pn']]){
          
        }; rm(sample.loop)
        
        
        ### Fix p-values etc.
        # Set up sig windows lists
        sentence.vs.list.sig.windows[[lock.loop]][[region.loop]] <- list()
        sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
        list.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
        
        # Reorganize p values
        sentence.vs.list[[region.loop]] <- 
          data.frame('p' = t(bind_rows(sentence.vs.list[[region.loop]])))
        sentence.vs.naming[[region.loop]] <- 
          data.frame('p' = t(bind_rows(sentence.vs.naming[[region.loop]])))
        list.vs.naming[[region.loop]] <- 
          data.frame('p' = t(bind_rows(list.vs.naming[[region.loop]])))
        
        ## Fill in sig windows
        for(alpha.loop in c(0.01, 0.05)){
          # alpha.loop = c(.01, .05)[1]
          
          alpha.label = paste0('alpha=',alpha.loop)
          
        ## Sentence vs. list
        sentence.vs.list[[region.loop]][,alpha.label] <- 
          as.numeric(sentence.vs.list[[region.loop]]$p < alpha.loop)
        sentence.vs.list.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
          get.significant.windows(sentence.vs.list[[region.loop]][,alpha.label],
                                  .sample.labels = row.names(sentence.vs.list[[region.loop]]),
                                  include.duration = TRUE,
                                  # .exclude.sig.durations.under.ms = 100,
                                  output.class = 'data.frame')
        
        ## Sentence vs. naming
        sentence.vs.naming[[region.loop]][,alpha.label] <- 
          as.numeric(sentence.vs.naming[[region.loop]]$p < alpha.loop)
        sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
          get.significant.windows(sentence.vs.naming[[region.loop]][,alpha.label],
                                  .sample.labels = row.names(sentence.vs.naming[[region.loop]]),
                                  include.duration = TRUE,
                                  # .exclude.sig.durations.under.ms = 100,
                                  output.class = 'data.frame')
        
        ## List vs. naming
        list.vs.naming[[region.loop]][,alpha.label] <- 
          as.numeric(list.vs.naming[[region.loop]]$p < alpha.loop)
        list.vs.naming.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
          get.significant.windows(list.vs.naming[[region.loop]][,alpha.label],
                                  .sample.labels = row.names(list.vs.naming[[region.loop]]),
                                  include.duration = TRUE,
                                  # .exclude.sig.durations.under.ms = 100,
                                  output.class = 'data.frame')
        
        }; rm(alpha.loop)
        
      }; rm(region.loop) # duration
      
      beep(2)
    }; rm(lock.loop)
    
    ## Save
    # Make a list of things you want to save:
    save.these <- c('region.means',
                    'region.ses',
                    'sample.labels',
                    'times',
                    'median.rt.samples',
                    'sentence.vs.list.sig.windows',
                    'sentence.vs.naming.sig.windows',
                    'list.vs.naming.sig.windows')
    
    # Save
    dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
    save(list = save.these,
         file = paste0(save.data.dir, 
                       'ecog wilcoxon task tests',
                       #'_',Sys.time(), # comment this line out if intentional save
                       '.RData'))
    
  } # if(file.exists(paste0(save.data.dir, 'ecog wilcoxon task tests.RData'))){
  
  
  
  ###
  ### Plots
  ###
  
  ### For region means and region ses, duplicate with sample labels
  region.means <- list('locked_to_production_onset' = region.means,
                       'locked_to_stimulus_onset' = region.means)
  region.ses <- list('locked_to_production_onset' = region.ses,
                     'locked_to_stimulus_onset' = region.ses)
  # Update stim-locked labels
  for(task.loop in tasks){
    # task.loop = tasks[1]
    for(region.loop in names(region.means[[1]][[task.loop]])){
      names(region.means$locked_to_stimulus_onset[[task.loop]][[region.loop]]) <-
        sample.labels$locked_to_stimulus_onset[[task.loop]]
      names(region.ses$locked_to_stimulus_onset[[task.loop]][[region.loop]]) <-
        sample.labels$locked_to_stimulus_onset[[task.loop]]
    }; rm(region.loop)
  }; rm(task.loop)
  
  
  
  ## Get y limits across prod-locked and stim-locked data
  y.limits <- list()
  regions.with.data <- names(which(sapply(region.means$locked_to_production_onset$pn, length) != 0))
  for(region.loop in regions.with.data){
    # region.loop = regions.with.data[1]
    y.limits[[region.loop]] <- c()
    
    for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
      # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[2]
      
      y.limits[[region.loop]] <-
        c(y.limits[[region.loop]],
          region.means[[lock.loop]]$pn[[region.loop]] + region.ses[[lock.loop]]$pn[[region.loop]],
          region.means[[lock.loop]]$pn[[region.loop]] - region.ses[[lock.loop]]$pn[[region.loop]],
          region.means[[lock.loop]]$lp[[region.loop]] + region.ses[[lock.loop]]$lp[[region.loop]],
          region.means[[lock.loop]]$lp[[region.loop]] - region.ses[[lock.loop]]$lp[[region.loop]],
          region.means[[lock.loop]]$sp[[region.loop]] + region.ses[[lock.loop]]$sp[[region.loop]],
          region.means[[lock.loop]]$sp[[region.loop]] - region.ses[[lock.loop]]$sp[[region.loop]])
    }; rm(lock.loop)
    
    y.limits[[region.loop]] <-
      c('min' = min(y.limits[[region.loop]]),
        'max' = max(y.limits[[region.loop]]))
  }; rm(region.loop)
  
  plot.regions <- c('pericentral','STG','IPL','MTG','IFG','MFG')
  plot.regions <- plot.regions[which(plot.regions %in% regions.with.data)]
  pub.regions.y.limits <- 
    c('min' = min(unlist(y.limits[plot.regions])),
      'max' = max(unlist(y.limits[plot.regions])))
  
  
  # Loop thru black and white backgrounds
  for(theme.loop in c('black','white')){
    # theme.loop = 'white'
    
    for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
      # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[1]
      
      # Loop through alpha thresholds for one-tailed significance
      for(alpha.label in names(sentence.vs.list.sig.windows[[1]][[1]])){
        # alpha.label = names(sentence.vs.list.sig.windows[[1]][[1]])[1]
        
        alpha.loop = as.numeric(gsub('alpha=','',alpha.label))
        
      ##
      ## Plot neural activity for each task by region
      ##
      
      regions.with.data <- names(which(sapply(region.means[[lock.loop]]$pn, length) != 0))
      
      ### All 3 comparisons plotted:
      
      # Save directory
      save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/by region/',alpha.label,'/3 comparisons/',theme.loop,'/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      for(region.loop in 1:length(regions.with.data)){
        # region.loop = 1
        current.region <- regions.with.data[region.loop]
        
        # Any electrodes in this region?
        if((length(region.means[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.means[[lock.loop]]$pn[[current.region]])))) &
           (length(region.ses[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.ses[[lock.loop]]$pn[[current.region]]))))){
          pdf(paste0(save.dir,current.region,'_',lock.loop,'_',alpha.label,'.pdf'),
              width = 7, height = 6.5)
          par(mfcol = c(1,1),
              oma = c(0,0,1,0),
              mar = c(0,0,3,0))
          plot.time.series(.y.values = list(region.means[[lock.loop]]$pn[[current.region]],
                                            region.means[[lock.loop]]$lp[[current.region]],
                                            region.means[[lock.loop]]$sp[[current.region]]),
                           .x.values = times[[lock.loop]][c('pn','lp','sp')],
                           .y.limits = y.limits[[current.region]],
                           .title = current.region,
                           .y.label = 'EGoG',
                           .x.label = 'time (ms)',
                           .colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                           .error.bars = list(region.ses[[lock.loop]]$pn[[current.region]],
                                              region.ses[[lock.loop]]$lp[[current.region]],
                                              region.ses[[lock.loop]]$sp[[current.region]]),
                           .sig.windows = list(
                             subset(list.vs.naming.sig.windows[[lock.loop]][[current.region]][[alpha.label]], duration > 100),
                             subset(sentence.vs.naming.sig.windows[[lock.loop]][[current.region]][[alpha.label]], duration > 100),
                             subset(sentence.vs.list.sig.windows[[lock.loop]][[current.region]][[alpha.label]], duration > 100)),
                           .sig.color = list(
                             colors[[theme.loop]]$task_line_plots[c('lp','pn'),'hex'],
                             colors[[theme.loop]]$task_line_plots[c('pn','sp'),'hex'],
                             colors[[theme.loop]]$task_line_plots[c('sp','lp'),'hex']),
                           .time.lock = lock.loop,
                           .zoom = 1.2,
                           .theme = theme.loop)
          add.text.line.multiple.colors(text.segments = c('picture naming','list production','sentence production'),
                                        text.colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                                        .outer = TRUE)
          dev.off()  
          
        } # if any elecs
      }; rm(region.loop)
      
      
      ### Just 2 comparisons plotted (sentences vs. {lists, words})
      
      # Save directory
      save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/by region/',alpha.label,'/2 comparisons/',theme.loop,'/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      for(region.loop in 1:length(regions.with.data)){
        # region.loop = 1
        current.region <- regions.with.data[region.loop]
        
        # Any electrodes in this region?
        if((length(region.means[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.means[[lock.loop]]$pn[[current.region]])))) &
           (length(region.ses[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.ses[[lock.loop]]$pn[[current.region]]))))){
          pdf(paste0(save.dir,current.region,'_',lock.loop,'_',alpha.label,'.pdf'),
              width = 7, height = 6.5)
          par(mfcol = c(1,1),
              oma = c(0,0,1,0),
              mar = c(0,0,3,0))
          plot.time.series(.y.values = list(region.means[[lock.loop]]$pn[[current.region]],
                                            region.means[[lock.loop]]$lp[[current.region]],
                                            region.means[[lock.loop]]$sp[[current.region]]),
                           .x.values = times[[lock.loop]][c('pn','lp','sp')],
                           .y.limits = y.limits[[current.region]],
                           .title = current.region,
                           .y.label = 'ECoG',
                           .x.label = 'time (ms)',
                           .colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                           .error.bars = list(region.ses[[lock.loop]]$pn[[current.region]],
                                              region.ses[[lock.loop]]$lp[[current.region]],
                                              region.ses[[lock.loop]]$sp[[current.region]]),
                           .sig.windows = list(
                             subset(sentence.vs.naming.sig.windows[[lock.loop]][[current.region]][[alpha.label]], duration > 100),
                             subset(sentence.vs.list.sig.windows[[lock.loop]][[current.region]][[alpha.label]], duration > 100)),
                           .sig.color = list(
                             colors[[theme.loop]]$task_line_plots['pn','hex'],
                             colors[[theme.loop]]$task_line_plots['lp','hex']),
                           .time.lock = lock.loop,
                           .zoom = 1.2,
                           .theme = theme.loop)
          add.text.line.multiple.colors(text.segments = c('picture naming','list production','sentence production'),
                                        text.colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                                        .outer = TRUE)
          dev.off()  
          
        } # if any elecs
      }; rm(region.loop)
      
      
      ##
      ## Just 1 comparison plotted (sentences vs. lists)
      ##
      
      # Save directory
      save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/by region/',alpha.label,'/1 comparison/',theme.loop,'/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      for(region.loop in 1:length(regions.with.data)){
        # region.loop = 1
        current.region <- regions.with.data[region.loop]
        
        # Any electrodes in this region?
        if((length(region.means[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.means[[lock.loop]]$pn[[current.region]])))) &
           (length(region.ses[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.ses[[lock.loop]]$pn[[current.region]]))))){
          pdf(paste0(save.dir,current.region,'_',lock.loop,'_',alpha.label,'.pdf'),
              width = 7, height = 6.5 * (3/4))
          # layout(matrix(c(1,2), nrow = 1, ncol = 2),
          #        widths = list('locked_to_production_onset' = c(3, 0),
          #                      'locked_to_stimulus_onset' = c(2.24,.76))[[lock.loop]])
          par(mfrow = c(1,1),
              oma = c(0,0,1,0),
              mar = c(0,0,3,0))
          plot.time.series(.y.values = list(region.means[[lock.loop]]$pn[[current.region]],
                                            region.means[[lock.loop]]$lp[[current.region]],
                                            region.means[[lock.loop]]$sp[[current.region]]),
                           .x.values = times[[lock.loop]][c('pn','lp','sp')],
                           #.y.limits = pub.regions.y.limits,
                           .y.limits.min.at.least = pub.regions.y.limits['min'],
                           .y.limits.max.at.least = pub.regions.y.limits['max'],
                           .title = current.region,
                           .y.label = 'ECoG',
                           .y.ticks = c(0,.5,1,1.5),
                           .x.ticks = list('locked_to_production_onset' = c(-1000, -500, 0, 500, 1000),
                                           'locked_to_stimulus_onset' = c(0, 500, 1000))[[lock.loop]],
                           .x.tick.labels = list('locked_to_production_onset' = c('-1000', '', '0', '', '1000'),
                                                 'locked_to_stimulus_onset' = c('0', '', '1000'))[[lock.loop]],
                           .x.label = 'time (ms)',
                           .x.limits = c(ifelse(lock.loop == 'locked_to_stimulus_onset', 0, -1000),
                                         ifelse(lock.loop == 'locked_to_stimulus_onset', 2000, 1000)),
                           show.t0 = ifelse(lock.loop == 'locked_to_stimulus_onset', FALSE, TRUE),
                           .colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                           .error.bars = list(region.ses[[lock.loop]]$pn[[current.region]],
                                              region.ses[[lock.loop]]$lp[[current.region]],
                                              region.ses[[lock.loop]]$sp[[current.region]]),
                           .sig.windows = list(
                             sentence.vs.list.sig.windows[[lock.loop]][[current.region]][[alpha.label]]),
                           .time.lock = lock.loop,
                           .sig.color = colors[[theme.loop]]$task_line_plots['sp','hex'],
                           .zoom = 1.2,
                           .theme = theme.loop)
          add.text.line.multiple.colors(text.segments = c('picture naming','list production','sentence production'),
                                        text.colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                                        .outer = TRUE)
          dev.off()  
          
        } # if any elecs
      }; rm(region.loop)
      
      
      ##
      ## Just 1 comparison plotted (sentences vs. lists) - PDF FOR PUBLICATION (no axis labels etc.)
      ##
      
      # Save directory
      save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/by region/',alpha.label,'/1 comparison - manuscript/',theme.loop,'/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      for(region.loop in 1:length(regions.with.data)){
        # region.loop = 1
        current.region <- regions.with.data[region.loop]
        
        # Any electrodes in this region?
        if((length(region.means[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.means[[lock.loop]]$pn[[current.region]])))) &
           (length(region.ses[[lock.loop]]$pn[[current.region]]) > 0) & 
           ((! any(is.na(region.ses[[lock.loop]]$pn[[current.region]]))))){
          
          current.y <- list(region.means[[lock.loop]]$pn[[current.region]],
                            region.means[[lock.loop]]$lp[[current.region]],
                            region.means[[lock.loop]]$sp[[current.region]])
          current.x <- times[[lock.loop]][c('pn','lp','sp')]
          current.error.bar <- list(region.ses[[lock.loop]]$pn[[current.region]],
                                    region.ses[[lock.loop]]$lp[[current.region]],
                                    region.ses[[lock.loop]]$sp[[current.region]])
          current.sig.window <- subset(sentence.vs.list.sig.windows[[lock.loop]][[current.region]][[alpha.label]],
                                       duration >= 100)
          if(lock.loop == 'locked_to_stimulus_onset'){
            max.plot.time.ms <- 1100
            current.y <- lapply(current.y, function(x){
              x[1:time.convert(100 + max.plot.time.ms, "times", "samples")]})
            current.x <- lapply(current.x, function(x){
              x[1:time.convert(100 + max.plot.time.ms, "times", "samples")]})
            current.error.bar <- lapply(current.error.bar, function(x){
              x[1:time.convert(100 + max.plot.time.ms, "times", "samples")]})
            current.sig.window$end.time <- sapply(current.sig.window$end.time, function(x){return(min(c(x, max.plot.time.ms)))})
            if(nrow(current.sig.window > 0)){
              current.sig.window$duration <- current.sig.window$end.time - current.sig.window$start.time  
            }
            current.sig.window <- subset(current.sig.window, duration >= 100)
          } # if(lock.loop == 'locked_to_stimulus_onset'){
          
          pdf(paste0(save.dir,current.region,'_',lock.loop,'_',alpha.label,'.pdf'),
              width = 8, height = 4)
          par(mfrow = c(1,1),
              oma = c(0,0,1,0),
              mar = c(0,0,3,0))
          
          plot.time.series(.y.values = current.y,
                           .x.values = current.x,
                           #.y.limits = pub.regions.y.limits,
                           .y.limits.min.at.least = pub.regions.y.limits['min'],
                           .y.limits.max.at.least = pub.regions.y.limits['max'],
                           # .title = current.region,
                           .y.label = '',#band.loop,
                           .y.ticks = list('locked_to_production_onset' = list('beta' = c(0,.5,1),
                                                                               'high_gamma' = c(0,1))[[band.loop]],
                                           'locked_to_stimulus_onset' = list('beta' = c(-.25, 0, .25, .5),
                                                                             'high_gamma' = c(0,1))[[band.loop]])[[lock.loop]],
                           .y.tick.labels = list('locked_to_production_onset' = list('beta' = c("0",'.5','1'),
                                                                                     'high_gamma' = c("0","1"))[[band.loop]],
                                                 'locked_to_stimulus_onset' = list('beta' = c('', '0', ' ', '.5'),
                                                                                   'high_gamma' = c("0","1"))[[band.loop]])[[lock.loop]],
                           .x.ticks = list('locked_to_production_onset' = c(-1000, -500, 0, 500, 1000),
                                           'locked_to_stimulus_onset' = c(0, 500, 1000))[[lock.loop]],
                           .x.tick.labels = list('locked_to_production_onset' = c('-1000', '', '0', '', '1000'),
                                                 'locked_to_stimulus_onset' = c('0', '', '1000'))[[lock.loop]],
                           .x.label = '',#time (ms)',
                           .x.limits = list('locked_to_production_onset' = c(-1600, 1100),
                                            'locked_to_stimulus_onset' = c(-100, 2600))[[lock.loop]],
                           show.t0 = TRUE,#(lock.loop == 'locked_to_production_onset'),
                           .horizontal.line.at = list('high_gamma' = NA,
                                                      'beta' = 0)[[band.loop]],
                           show.y.axis = (lock.loop == 'locked_to_stimulus_onset'),
                           .colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                           .error.bars = current.error.bar,
                           .sig.windows = list(current.sig.window),
                           .sig.color = colors[[theme.loop]]$task_line_plots['sp','hex'],
                           .zoom = 1.8,
                           .margin = c(3,3,0,.5),
                           .theme = theme.loop,
                           .background = ifelse(theme.loop == 'white', rgb(1,1,1,0), theme.loop))
          dev.off()  
          
        } # if any elecs
      }; rm(region.loop)
      
    } # alpha.label and alpha.loop
      
      ##
      ## Plot ECoG for each region by task
      ## 
      
      
      # Save directory
      save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/electrodes/ecog by task/mean/',theme.loop,'/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      for(task.loop in tasks){
        # task.loop = tasks[1]
        
        pdf(paste0(save.dir,task.loop,'_',lock.loop,'.pdf'),
             width = 7, height = 6.5)
        par(mfcol = c(1,1),
            oma = c(0,0,1,0),
            mar = c(0,0,3,0))
        plot.time.series(.y.values = region.means[[lock.loop]][[task.loop]][plot.regions],
                         .x.values = times[[lock.loop]][task.loop],
                         .title = task.loop,
                         .colors = colors[[theme.loop]]$rois[plot.regions,'hex'],
                         #.y.limits = c(0, 1.6),
                         .y.label = 'ECoG',
                         .x.label = 'time (ms)',
                         .error.bars = region.ses[[lock.loop]][[task.loop]][plot.regions],
                         .time.lock = lock.loop,
                         .zoom = 1.2,
                         .theme = theme.loop)
        
        add.text.line.multiple.colors(text.segments = plot.regions,
                                      text.colors = colors[[theme.loop]]$rois[plot.regions,'hex'],
                                      .outer = TRUE)
        dev.off()  
      }
      
    }# lock.loop
    
    
    
    ##
    ## Just 1 comparison plotted (sentences vs. lists) - split by median RT
    ##
    
    # Loop through alpha thresholds for one-tailed significance
    for(alpha.label in names(sentence.vs.list.sig.windows[[1]][[1]])){
      # alpha.label = names(sentence.vs.list.sig.windows[[1]][[1]])[1]
      
      alpha.loop = as.numeric(gsub('alpha=','',alpha.label))
      
    
    # Save directory
    save.dir <- paste0(path,'analysis/R/task ECoG comparisons/ROI wilcox tests and squiggle plots/output/1b - warped data/figures/',band.loop,'/by region/',alpha.label,'/1 comparison - split by median RT/',theme.loop,'/')
    dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
    
    for(lock.loop in c('locked_to_stimulus_onset', 'locked_to_production_onset')){
      # lock.loop = c('locked_to_stimulus_onset', 'locked_to_production_onset')[1]
      
      regions.with.data <- names(which(sapply(region.means[[lock.loop]]$pn, length) != 0))
      
    for(region.loop in 1:length(regions.with.data)){
      # region.loop = 1
      current.region <- regions.with.data[region.loop]
      
      # Any electrodes in this region?
      if((length(region.means$locked_to_production_onset$pn[[current.region]]) > 0) & 
         ((! any(is.na(region.means$locked_to_production_onset$pn[[current.region]])))) &
         (length(region.ses$locked_to_production_onset$pn[[current.region]]) > 0) & 
         ((! any(is.na(region.ses$locked_to_production_onset$pn[[current.region]]))))){
        
          
        # Get y-vals
        current.ys <- list()
        current.error.bars <- list()
        for(task.loop in c('pn','lp','sp')){
          # task.loop = c('pn','lp','sp')[1]
          
          # Y-values
          current.ys[[task.loop]] <-
            region.means[[lock.loop]][[task.loop]][[current.region]]
          # Error bars
          current.error.bars[[task.loop]] <-
            region.ses[[lock.loop]][[task.loop]][[current.region]]
          
          # Subset
          if(lock.loop == 'locked_to_stimulus_onset'){
            ## If locked to stimulus onset
            # Ys
            current.ys[[task.loop]] <- current.ys[[task.loop]][
              which(time.convert(names(current.ys[[task.loop]]), "sample.labels", "samples") <
                      (floor(median.rt.samples[task.loop]/2)))]
            # Error bars
            current.error.bars[[task.loop]] <- current.error.bars[[task.loop]][
              which(time.convert(names(current.error.bars[[task.loop]]), "sample.labels", "samples") <
                      (floor(median.rt.samples[task.loop]/2)))]  
          }else{
            ## If locked to production onset
            # Ys
            current.ys[[task.loop]] <- current.ys[[task.loop]][
              which(time.convert(names(current.ys[[task.loop]]), "sample.labels", "samples") >
                      (-floor(median.rt.samples[task.loop]/2)))]
            # Error bars
            current.error.bars[[task.loop]] <- current.error.bars[[task.loop]][
              which(time.convert(names(current.error.bars[[task.loop]]), "sample.labels", "samples") >
                      (-floor(median.rt.samples[task.loop]/2)))]  
          }
          
        }; rm(task.loop)
        
        ## Sig windows
        # Get sig windows just for first half (stim-locked) and second half (prod-locked)
        current.sig.windows <- 
          sentence.vs.list.sig.windows[[lock.loop]][[current.region]][[alpha.label]]
        if(lock.loop == 'locked_to_stimulus_onset'){
          ## If locked to stimulus onset
        current.sig.windows$end.time <- 
          sapply(current.sig.windows$end.time, function(x){
            min(x, time.convert(floor(median.rt.samples['lp'] / 2), "samples", "times"))})
        }else{
          ## If locked to production onset
          current.sig.windows$start.time <- 
            sapply(current.sig.windows$start.time, function(x){
              max(x, time.convert(-ceiling(median.rt.samples['lp'] / 2), "samples", "times"))})
        }
        
        # Update durations
        if(nrow(current.sig.windows) > 0){
          current.sig.windows$duration <- current.sig.windows$end.time - current.sig.windows$start.time
          current.sig.windows <- subset(current.sig.windows, duration > 0)
        }
        
        
        pdf(paste0(save.dir,current.region,'_',lock.loop,'_',alpha.label,'.pdf'),
            width = 7, height = 4)
        par(mfrow = c(1,1),
            oma = c(0,0,0,0))
        plot.time.series(.y.values = current.ys,
                         .x.values = lapply(current.ys, names),
                         #.y.limits = pub.regions.y.limits,
                         .y.limits.min.at.least = pub.regions.y.limits['min'],
                         .y.limits.max.at.least = pub.regions.y.limits['max'],
                         # .title = current.region,
                         .y.label = '',#band.loop,
                         .y.ticks = c(0,.5,1),
                         .x.ticks = list('locked_to_production_onset' = c(-500, 0, 500, 1000),
                                         'locked_to_stimulus_onset' = c(0, 500))[[lock.loop]],
                         .x.label = '',#time (ms)',
                         .x.limits = list('locked_to_production_onset' = c(-600, 1100),
                                          'locked_to_stimulus_onset' = c(-100, 1600))[[lock.loop]],
                         show.t0 = TRUE,#(lock.loop == 'locked_to_production_onset'),
                         show.y.axis = (lock.loop == 'locked_to_stimulus_onset'),
                         .colors = colors[[theme.loop]]$task_line_plots[names(current.ys),'hex'],
                         .error.bars = current.error.bars,
                         .sig.windows = list(subset(current.sig.windows, duration > 100)),
                         .sig.color = colors[[theme.loop]]$task_line_plots['sp','hex'],
                         .zoom = 1.8,
                         .margin = c(3,3,0,0),
                         .theme = theme.loop,
                         .background = ifelse(theme.loop == 'white', rgb(1,1,1,0), theme.loop))
        # add.text.line.multiple.colors(text.segments = c('picture naming','list production','sentence production'),
        #                               text.colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
        #                               .outer = TRUE)
        dev.off()
        
        
      } # if any elecs
    }; rm(region.loop)
      
    }; rm(lock.loop)
    
  } # alpha.label and alpha.loop
    
    
    
    
  }# theme.loop
  
 
  
} # band.loop

#stopCluster(cl)

message('Script completed successfully. ',Sys.time())
