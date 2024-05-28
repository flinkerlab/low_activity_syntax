### Correlate RSIs and high gamma across electrodes at time when RSI peaks
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
# Adjust transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Get min and max
source(paste0(path,'/analysis/R/functions/min_max.R'))


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


### Clean up
keep <- c(ls(), 
          'keep.all.this', 
          'ecog.freq.loop', 
          'rsa.freq.loop',
          'alpha.hps',
          'sig.hp.loop')

### Loop thru stuff
for(ecog.freq.loop in c('high_gamma','beta')){ ## UNCOMMENT
  # ecog.freq.loop = c('high_gamma','beta')[1] ## UNCOMMENT
  
  for(rsa.freq.loop in c('high_gamma','beta')){ ## UNCOMMENT
    # rsa.freq.loop = c('high_gamma','beta')[1] ## UNCOMMENT
    
    
    # Define alpha hyperparameters to loop thru
    alpha.hps <- c('alpha=.05_min=100ms',
                   'alpha=.01_min=50ms')
  
    # Clean up
    rm(list = ls()[which(!ls() %in% keep)])
    gc()
    
    # Directories
    output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/', ecog.freq.loop,'/')
    ecog.input.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',ecog.freq.loop, '/')
    rsa.input.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',rsa.freq.loop, '/')
    
    
    
    ##
    ## Load all data
    ##
    
    ### Get median.rt.samples
    load(paste0(path, 'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
    
    
    ### Elec means - high gamma
    load(paste0(path,'analysis/R/electrode means/means - warped data - multiband/output - ',ecog.freq.loop,'/data/warped electrode data means, SDs, and SEs for all 3 tasks.RData')) # loads elec.means, elec.sds, and elec.ses
    
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
    elec.sds <- lapply(elec.sds, bind_cols)
    elec.ses <- lapply(elec.ses, bind_cols)
    
    
    ### Load unadjusted smoothed z-values
    unadjusted.ts.path <- paste0(rsa.input.path,'data/1c - stats - shuffled RSA/real data stats/')
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
    
    ### Load adjusted smoothed z-values
    adjusted.ts.path <- paste0(rsa.input.path,'data/1c - stats - shuffled RSA/real data stats/')
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
    
    
    ### Loop through significance thresholds for elecs
    keep <- c(ls(), "keep", "sig.hp.loop")
    # for(sig.hp.loop in alpha.hps){ ## UNCOMMENT
      sig.hp.loop = alpha.hps[1] ## UNCOMMENT
        
        # Clean up
        rm(list = ls()[which(!ls() %in% keep)])
        gc()
        
        current.sig.elecs <- zs.sig.elecs[[sig.hp.loop]]
        
        ### Get HGP at time of max adjusted t for each sig elec
        ecog.at.max.adjusted.z <- list()
        for(term.loop in names(current.sig.elecs)){
          # term.loop = names(current.sig.elecs)[1]
          
          # Set up dataframe for this term to store max t-adjusted and HGA for each sig elec for this term
          ecog.at.max.adjusted.z[[term.loop]] <- data.frame('elec' = current.sig.elecs[[term.loop]],
                                                            'max.t.adjusted' = NA,
                                                            'hga' = NA)
          rownames(ecog.at.max.adjusted.z[[term.loop]]) <- current.sig.elecs[[term.loop]]
          
          # Add colors for whether or not elec is in sensorimotor cortex
          ecog.at.max.adjusted.z[[term.loop]]$sensorimotor.coloring <- rgb(.75,.75,.75)
          ecog.at.max.adjusted.z[[term.loop]]$region_clinical <- 
            elec.info[ecog.at.max.adjusted.z[[term.loop]]$elec, 'region_clinical']
          ecog.at.max.adjusted.z[[term.loop]]$sensorimotor.coloring[ecog.at.max.adjusted.z[[term.loop]]$region_clinical %in% c('precentral', 'postcentral')] <-
            colors$black$rainbow_bright[c('green'),'hex']
          ecog.at.max.adjusted.z[[term.loop]]$sensorimotor.coloring[ecog.at.max.adjusted.z[[term.loop]]$region_clinical %in% c('cSTG','ctx-lh-superiortemporal','mSTG','rSTG')] <-
            colors$black$rainbow_bright[c('pink'),'hex']
          
          
          # Loop thru elecs
          for(sig.elec.loop in rownames(ecog.at.max.adjusted.z[[term.loop]])){
            # sig.elec.loop = rownames(ecog.at.max.adjusted.z[[term.loop]])[1]
            
            # ## Look only at elecs sig in planning window? (Doesn't change much)
            # current.sig.windows <- zs.sig.windows[[sig.elec.loop]][[term.loop]][[sig.hp.loop]]
            # # Set latest sig time to 0
            # current.sig.windows$end.time <- sapply(current.sig.windows$end.time, function(x){min(x, 0)})
            # current.sig.windows$duration <- current.sig.windows$end.time - current.sig.windows$start.time + 1
            # current.sig.windows <- current.sig.windows[current.sig.windows$duration >= 100,]
            # 
            # # If any sig windows
            # if(nrow(current.sig.windows) > 0){
            
              ## Get the sample.label where t.adjusted peaks
              max.index <- which.max(zs.adjusted[[sig.elec.loop]][,term.loop])[1]
              max.sample.label <- rownames(zs.adjusted[[sig.elec.loop]])[max.index]
              
              # Store t.adjusted and HGA
              ecog.at.max.adjusted.z[[term.loop]][sig.elec.loop,'max.t.adjusted'] <- 
                zs.adjusted[[sig.elec.loop]][max.sample.label, term.loop]
              ecog.at.max.adjusted.z[[term.loop]][sig.elec.loop,'hga'] <-
                elec.means[['sp']][max.sample.label, sig.elec.loop]
            
            # }else{ # if(nrow(current.sig.windows) > 0){
            #   
            #   # Remove elec from dataframe
            #   ecog.at.max.adjusted.z[[term.loop]] <- 
            #     ecog.at.max.adjusted.z[[term.loop]][-which(rownames(ecog.at.max.adjusted.z[[term.loop]]) == sig.elec.loop),]
            #   
            # } # if(nrow(current.sig.windows) > 0){}else{}
            
          }; rm(sig.elec.loop)
          
        }#; rm(term.loop)
        
        ### Remove outliers 
        # Peak < p = .05
        for(term.loop in names(current.sig.elecs)){
          # term.loop = names(current.sig.elecs)[1]
          
          ecog.at.max.adjusted.z[[term.loop]] <- 
            ecog.at.max.adjusted.z[[term.loop]][! ecog.at.max.adjusted.z[[term.loop]]$max.t.adjusted < qnorm(.95) ,]
        }#; rm(term.loop)
        
        # Using boxplot()'s criterion (I think 1.5 * IQI + {3rd Q})
        
          for(term.loop in names(current.sig.elecs)){
            # term.loop = names(current.sig.elecs)[1]
            
            ecog.at.max.adjusted.z[[term.loop]] <- 
              ecog.at.max.adjusted.z[[term.loop]][! ecog.at.max.adjusted.z[[term.loop]]$hga %in% 
                                                    boxplot(ecog.at.max.adjusted.z[[term.loop]]$hga, 
                                                            plot = FALSE)$out,]
            ecog.at.max.adjusted.z[[term.loop]] <- 
              ecog.at.max.adjusted.z[[term.loop]][! ecog.at.max.adjusted.z[[term.loop]]$max.t.adjusted %in% 
                                                    boxplot(ecog.at.max.adjusted.z[[term.loop]]$max.t.adjusted, 
                                                            plot = FALSE)$out,]
          }#; rm(term.loop)
        
        
        
        # Range of all values
        range <- bind_rows(ecog.at.max.adjusted.z)[,c('max.t.adjusted','hga')]
        range <- lapply(range, function(x){c(min(x), max(x))})
        
        
        ### Plot!
        # for(theme.loop in c('black','white')){ ## UNCOMMENT
          theme.loop = c('black','white')[2] ## UNCOMMENT
          
          # Plot parameters
          term.colors <- colors[[theme.loop]]$rsa_term_line_plots$hex
          names(term.colors) <- rownames(colors[[theme.loop]]$rsa_term_line_plots)
          text.size.big <- 1.8
          text.size.med <- 1.6
          text.size.small <- 1.4
          zoom <- 1.8
          
          # Plot.dir
          plot.dir <- 
            paste0(output.path, 'figures/1k - ECoG RSA relationship/1k_a - t-value scatterplots/main RSA model/',
                   theme.loop,
                   '/individual axis limits/',
                   sig.hp.loop,'/',
                   '/rsa terms - ',rsa.freq.loop,'/')
          plot.dir.same.axes <- 
            paste0(output.path, 
                   'figures/1k - ECoG RSA relationship/1k_a - t-value scatterplots/main RSA model/',
                   theme.loop,
                   '/same axis limits/',
                   sig.hp.loop,'/',
                   '/rsa terms - ',rsa.freq.loop,'/',
                   '/term colored/')
          plot.dir.region.colored <- 
            paste0(output.path, 
                   'figures/1k - ECoG RSA relationship/1k_a - t-value scatterplots/main RSA model/',
                   theme.loop,
                   '/same axis limits/',
                   sig.hp.loop,'/',
                   '/rsa terms - ',rsa.freq.loop,'/',
                   '/region colored/')
          plot.dir.same.axes.square <- 
            paste0(output.path, 
                   'figures/1k - ECoG RSA relationship/1k_a - t-value scatterplots/main RSA model/',
                   theme.loop,
                   '/same axis limits/',
                   sig.hp.loop,'/',
                   '/rsa terms - ',rsa.freq.loop,'/',
                   '/term colored - square/')
          dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(plot.dir.same.axes, showWarnings = FALSE, recursive = TRUE)
          dir.create(plot.dir.same.axes.square, showWarnings = FALSE, recursive = TRUE)
          dir.create(plot.dir.region.colored, showWarnings = FALSE, recursive = TRUE)
          
          for(term.loop in terms){
            # term.loop = terms[2]
            
            current.data <- ecog.at.max.adjusted.z[[term.loop]]
            
            for(plot.loop in c('same.range_term.colored',
                               'same.range_region.colored',
                               'same.range_term.colored.square',
                               'different.range_term.colored')){
              # plot.loop = c('same.range_term.colored','same.range_region.colored','different.range_term.colored')[1]
              
              # Current model 
              model <- with(current.data, lm(hga ~ max.t.adjusted))
              current.p <- summary(model)$coefficients[2,4]
              
              if(plot.loop == 'same.range_term.colored'){
                ## Term colored plot
                pdf(paste0(plot.dir.same.axes, term.loop,' - ',rsa.freq.loop,' - p=',round(current.p, 3),'.pdf'),
                    width = 4, height = 4)
                par(bg = ifelse(theme.loop == 'black', 'black', rgb(1,1,1,0)),
                    pty = 's',
                    # mar = c(6,6,3,1) * zoom,
                    mar = c(3,3,0,.5) * zoom,
                    oma = c(0,0,1,0))
                plot(x = current.data$max.t.adjusted,
                     y = current.data$hga,
                     xlab = '',
                     xlim = range[['max.t.adjusted']],#c(range[['max.t.adjusted']][1] * .95, range[['max.t.adjusted']][2]),
                     ylim = range[['hga']],
                     ylab = '',
                     pch = 20,
                     cex = 1.6 * zoom,
                     col = adjust.transparency(term.colors[term.loop], .5),
                     bty = 'n',
                     xaxt = 'n',
                     yaxt = 'n')
              } # if(plot.loop == same.range_term.colored){
              if(plot.loop == 'same.range_term.colored.square'){
                ## Term colored plot
                pdf(paste0(plot.dir.same.axes.square, term.loop,' - ',rsa.freq.loop,' - p=',round(current.p, 3),'.pdf'),
                    width = 5, height = 5)
                par(bg = ifelse(theme.loop == 'black', 'black', rgb(1,1,1,0)),
                    pty = 's',
                    # mar = c(6,6,3,1) * zoom,
                    mar = c(3,3,0,.5) * zoom,
                    oma = c(0,0,1,0))
                plot(x = current.data$max.t.adjusted,
                     y = current.data$hga,
                     xlab = '',
                     xlim = range[['max.t.adjusted']],
                     ylim = range[['hga']],
                     ylab = '',
                     pch = 20,
                     cex = 1.6 * zoom,
                     col = adjust.transparency(term.colors[term.loop], .5),
                     bty = 'n',
                     xaxt = 'n',
                     yaxt = 'n')
              } # if(plot.loop == same.range_term.colored){
              if(plot.loop == 'same.range_region.colored'){
                ## Region colored plot
                pdf(paste0(plot.dir.region.colored, term.loop,' - ',rsa.freq.loop,' - p=',round(current.p, 3),'.pdf'),
                    width = 6.5, height = 4)
                par(bg = ifelse(theme.loop == 'black', 'black', rgb(1,1,1,0)),
                    pty = 's',
                    # mar = c(6,6,3,1) * zoom,
                    mar = c(3,3,0,.5) * zoom,
                    oma = c(0,0,1,0))
                plot(x = current.data$max.t.adjusted,
                     y = current.data$hga,
                     xlab = '',
                     xlim = range[['max.t.adjusted']],
                     ylim = range[['hga']],
                     ylab = '',
                     pch = 20,
                     cex = 1.6 * zoom,
                     col = adjust.transparency(current.data$sensorimotor.coloring, .5),
                     bty = 'n',
                     xaxt = 'n',
                     yaxt = 'n')
              } # if(plot.loop == 'same.range_region.colored'){
              if(plot.loop == 'different.range_term.colored'){
                pdf(paste0(plot.dir, term.loop,' - ',rsa.freq.loop,' - p=',round(current.p, 3),'.pdf'),
                    width = 6.5, height = 4)
                par(bg = ifelse(theme.loop == 'black', 'black', rgb(1,1,1,0)),
                    pty = 's',
                    # mar = c(6,6,3,1) * zoom,
                    mar = c(3,3,0,.5) * zoom,
                    oma = c(0,0,1,0))
                plot(x = current.data$max.t.adjusted,
                     y = current.data$hga,
                     xlab = '',
                     ylab = '',
                     pch = 20,
                     cex = 1.6 * zoom,
                     col = adjust.transparency(term.colors[term.loop], .5),
                     bty = 'n',
                     xaxt = 'n',
                     yaxt = 'n')
              } # if(plot.loop == 'different.range_term.colored'){
              abline(model,
                     lwd = 2 * zoom,
                     col = ifelse(theme.loop == 'black', 'white', 'black'))
              
              
              
              ### Axes
              if(plot.loop == 'different.range_term.colored'){
                x.range <- c(ceiling(min(current.data$max.t.adjusted)),
                             floor(max(current.data$max.t.adjusted)))
                x.at <- x.range[1]:x.range[2]
                x.at.labels <- as.character(x.at)
                x.at.labels[! x.at %in% c(x.at[1], x.at[length(x.at)])] <- ""
                y.range <- c(ceiling(min(current.data$hga)),
                             floor(max(current.data$hga)))
                if(y.range[1] == y.range[2]){
                  y.range <- c(y.range[1] - .5,
                               y.range[2] + .5)
                }
                y.at <- y.range[1]:y.range[2]
                y.at.labels <- as.character(y.at)
                y.at.labels[! y.at %in% c(y.at[1], y.at[length(y.at)])] <- ""
                axis(side = 1,
                     at = x.range,
                     labels = x.range,
                     las = 1,
                     tck = -.025 * zoom, # length of tick
                     padj = .6 * zoom, # distance between tick and label
                     lwd = 1.5 * zoom,
                     lwd.ticks = 1.5 * zoom,
                     cex.axis = text.size.med * zoom,
                     col = ifelse(theme.loop == 'black', 'white', 'black'),
                     col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
                axis(side = 2,
                     at = y.range,
                     labels = y.range,
                     las = 0,
                     tck = -.025 * zoom, # length of tick
                     padj = -.45 * zoom, # distance between tick and label
                     lwd = 1.5 * zoom,
                     lwd.ticks = 1.5 * zoom,
                     cex.axis = text.size.med * zoom,
                     col = ifelse(theme.loop == 'black', 'white', 'black'),
                     col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
              }else{ # if(plot.loop == 'different.range_term.colored'){
                
                ## Axes
                x.range <- c(ceiling(range$max.t.adjusted[1]),
                             floor(range$max.t.adjusted[2]))
                y.range <- c(ceiling(range$hga[1]),
                             floor(range$hga[2]))
                if(y.range[1] == y.range[2]){
                  y.range <- c(y.range[1] - .5,
                               y.range[2] + .5)
                }
                axis(side = 1,
                     at = x.range,
                     labels = x.range,
                     las = 1,
                     tck = -.025 * zoom, # length of tick
                     padj = .6 * zoom, # distance between tick and label
                     lwd = 1.5 * zoom,
                     lwd.ticks = 1.5 * zoom,
                     cex.axis = text.size.med * zoom,
                     col = ifelse(theme.loop == 'black', 'white', 'black'),
                     col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
                axis(side = 2,
                     at = y.range,
                     labels = y.range,
                     las = 0,
                     tck = -.025 * zoom, # length of tick
                     padj = -.45 * zoom, # distance between tick and label
                     lwd = 1.5 * zoom,
                     lwd.ticks = 1.5 * zoom,
                     cex.axis = text.size.med * zoom,
                     col = ifelse(theme.loop == 'black', 'white', 'black'),
                     col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
              }
              dev.off()
              
            }#; rm(plot.loop)
          }#; rm(term.loop)
        }#; rm(theme.loop)
        
    }#; rm(sig.hp.loop)
    
  }#; rm(rsa.freq.loop)
}#; rm(ecog.freq.loop)

# Finish!
message('Script completed successfully. ',Sys.time())







