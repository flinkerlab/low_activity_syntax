### Compare RSA outcomes on real data to shuffled data to derive RSIs/z-scored t-values for each representation type
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

# Clean up
rm(list=ls())
cat("\014")
message("Begin stats on distributions of shuffled RSA results. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 4
  n.cores.to.use.less = 2
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 6
    n.cores.to.use.less = 3
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 12
    n.cores.to.use.less = 8
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 44
  n.cores.to.use.less = 6
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


### Former loops
band.loop = c('high_gamma','beta')[2]
include.hga.loop = c(TRUE,FALSE)[1]

# Directories
read.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/', band.loop, '/')
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/', band.loop, '/')



##
## Load shuffle data
##

shuffle.ts.path <- paste0(output.path, '/data/1b - elec ts - shuffled trial data/shuffle ts smooth/')
shuffle.ts.files <- list.files(shuffle.ts.path)


### Read data
## Set up parallel processing
# Close any old parallel backends
unregister_dopar()
# Set up parallel workers in case caret wants to use it
cl <- makeCluster(n.cores.to.use, type = "FORK")
registerDoParallel(cl)

shuffle.ts.smooth <- 
  foreach(shuffle.loop = shuffle.ts.files) %dopar% {
    # shuffle.loop = shuffle.ts.files[1]
    
    load(paste0(shuffle.ts.path, shuffle.loop))
    return(sp.shuffle.ts.smooth)
    
  } # shuffle.loop

## Close parallel backend
stopCluster(cl)
unregister_dopar()

message("Done loading ",length(shuffle.ts.smooth)," shuffle stats files! ",Sys.time())
if(length(shuffle.ts.smooth) != 1000){message("WARNING!! Number of files loaded is ",length(shuffle.ts.smooth),", not 1000 as expected.")}
beep()


### Get summary stats
start.time <- Sys.time()
message("Getting summary stats of shuffled data. ", start.time)

# Storage
shuffle.means <- list()
shuffle.sds <- list()

# Get summary stats
for(elec.loop in names(shuffle.ts.smooth[[1]])){
  # elec.loop = names(shuffle.ts.smooth[[1]])[1]
    
    # Get data
    current.data <- lapply(shuffle.ts.smooth, function(x){x[[elec.loop]]})
    
    # Get means and SDs across loops
    shuffle.means[[elec.loop]] <- elementwise.matrix.apply(current.data, .function = "mean")
    shuffle.sds[[elec.loop]] <- elementwise.matrix.apply(current.data, .function = "sd")
    
}; rm(elec.loop)
message("Shuffled data summary stats done.  Duration: ",round(Sys.time() - start.time, 2)); rm(start.time)
beep()


### Clean up
rm(shuffle.ts.smooth)




##
## Load real smoothed t-values
##

smoothed.ts.path <- paste0(read.path,'data/1a - elec ts/')
smoothed.ts.file <- 'elec RSA multiple regression t-values for PN and SP.RData'
load(paste0(smoothed.ts.path, smoothed.ts.file))


### Get stats
# Significance
significance.hps <- list(
  'alpha=.05_min=100ms' = c('sig.threshold' = qnorm(1 - .05),
                            'min.sig.window' = 100),
  'alpha=.01_min=50ms' = c('sig.threshold' = qnorm(1 - .01),
                            'min.sig.window' = 50),
  'alpha=.01_min=100ms' = c('sig.threshold' = qnorm(1 - .01),
                            'min.sig.window' = 100),
  'alpha=.001_min=50ms' = c('sig.threshold' = qnorm(1 - .001),
                            'min.sig.window' = 50))

# Storage
ts.z <- list()
zs.sig.windows <- list()

# Loop thru elecs getting significant windows
for(elec.loop in names(shuffle.means)){
  # elec.loop = names(shuffle.means)[1]
  
  # Storage
  zs.sig.windows[[elec.loop]] <- list()
  ts.z[[elec.loop]] <- 
    data.frame(apply(sp.ts.smooth[[elec.loop]][,sp.model.terms[sp.model.terms != "diff.log.rt"]], 
                     c(1,2), function(x){return(NA)}))
  
  for(term.loop in colnames(ts.z[[elec.loop]])){
    # term.loop = sp.model.terms[1]
    
    # Get z-score
    ts.z[[elec.loop]][,term.loop] <-
      (sp.ts.smooth[[elec.loop]][,term.loop] - 
         shuffle.means[[elec.loop]][,term.loop]) /
      shuffle.sds[[elec.loop]][,term.loop]
    
    # Get sig windows at alpha = .05 for >= 100ms
    zs.sig.windows[[elec.loop]][[term.loop]] <- list()
    for(significance.loop in names(significance.hps)){
      zs.sig.windows[[elec.loop]][[term.loop]][[significance.loop]] <- 
        get.significant.windows(as.numeric(ts.z[[elec.loop]][,term.loop] > significance.hps[[significance.loop]]['sig.threshold']),
                                .sample.labels = rownames(ts.z[[elec.loop]]),
                                output.class = 'data.frame',
                                include.duration = TRUE,
                                .exclude.times.before.ms = time.convert(-median.rt.samples['sp'], "sample.labels", "times"),
                                .exclude.times.after.ms = 500,
                                .exclude.sig.durations.under.ms = significance.hps[[significance.loop]]['min.sig.window'])
    }; rm(significance.loop)
    
  }; rm(term.loop)
}; rm(elec.loop)


### Get sig elecs
zs.sig.elecs <- list()
n.sig.elecs <- list()
for(significance.loop in names(significance.hps)){
  zs.sig.elecs[[significance.loop]] <- list()
  
  for(term.loop in names(zs.sig.windows[[1]])){
    # term.loop = names(zs.sig.windows[[1]])[1]
    
    zs.sig.elecs[[significance.loop]][[term.loop]] <-
      sapply(zs.sig.windows, function(x){
        nrow(x[[term.loop]][[significance.loop]]) > 0
      })
    zs.sig.elecs[[significance.loop]][[term.loop]] <- 
      names(zs.sig.elecs[[significance.loop]][[term.loop]][zs.sig.elecs[[significance.loop]][[term.loop]]])
  }; rm(term.loop)
  
  ## Number of sig elecs
  n.sig.elecs[[significance.loop]] <- data.frame(t(sapply(zs.sig.elecs[[significance.loop]], length)))
  
}; rm(significance.loop)


###
### Adjust z-scored t-values
###

# Define function for squishing data: reduce extreme low and high vals
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
  out[out < 1] <- exp(out[out < 1]) / exp(1) # Reduce non-significant values
  out[out < 1] <- out[out < 1] ^ 3 # Even further reduce non-significant values
  out <- out * scaling.factor # Unscale  
  out[na.indices] <- NA
  return(out)
}

# Adjust
zs.adjusted <- list()
for(elec.loop in names(ts.z)){
  # elec.loop = names(ts.z)[1]
  
  # Storage
  zs.adjusted[[elec.loop]] <- 
    data.frame(apply(ts.z[[elec.loop]], 2, squish))
  
}; rm(elec.loop)



###
### Save
###


### Shuffle dist
save.these.shuffle.summary.stats <- c(
  'shuffle.means',
  'shuffle.sds'
)
save.shuffle.summary.stats.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/shuffle distribution summary stats/')
dir.create(save.shuffle.summary.stats.path, showWarnings = FALSE, recursive = TRUE)
save(list = save.these.shuffle.summary.stats,
     file = paste0(save.shuffle.summary.stats.path, 'summary stats for RSA shuffle distributions and z-scored ts.RData'))


### Z-scores, adjusted z's, and stats
# Z-scores
save.stats.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/real data stats/')
dir.create(save.stats.path, showWarnings = FALSE, recursive = TRUE)
save(ts.z,
     file = paste0(save.stats.path, 'just z-scores.RData'))

# Everything else
save.these <- c(
  'zs.adjusted',
  'zs.sig.windows',
  'zs.sig.elecs',
  'significance.hps'
)
save(list = save.these,
     file = paste0(save.stats.path, 'adjusted zs with sig windows and elecs.RData'))


### Number of sig elecs
save.n.sig.elecs.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/number of sig elecs/')
dir.create(save.n.sig.elecs.path, showWarnings = FALSE, recursive = TRUE)
save(n.sig.elecs,
     file = paste0(save.n.sig.elecs.path, 'n sig elecs - real data.RData'))



### Clean up
rm(ts.z, zs.adjusted, zs.sig.windows, zs.sig.elecs, sp.ts.smooth)



###
### Translate each shuffle loop's data into z-scores too to plug into NMF
###


# Close any old parallel backends
unregister_dopar()
# Set up parallel workers in case caret wants to use it
cl <- makeCluster(n.cores.to.use.less, type = "FORK")
registerDoParallel(cl)

temp.n.sig.elecs.in.shuffle.data <-
  foreach(shuffle.loop = 1:length(shuffle.ts.files)) %dopar% {
    # shuffle.loop = 1
    
    # Current shuffle data
    attach(paste0(shuffle.ts.path, shuffle.ts.files[shuffle.loop]))
    current.shuffle.ts <- sp.shuffle.ts.smooth
    detach()
    
    # Storage
    shuffle.ts.z <- list()
    shuffle.zs.sig.windows <- list()
    
    # Loop thru elecs getting significant windows
    for(elec.loop in names(shuffle.means)){
      # elec.loop = names(shuffle.means)[1]
      
      # Storage
      shuffle.zs.sig.windows[[elec.loop]] <- list()
      shuffle.ts.z[[elec.loop]] <- 
        data.frame(apply(current.shuffle.ts[[elec.loop]][,sp.model.terms[sp.model.terms != "diff.log.rt"]], 
                         c(1,2), function(x){return(NA)}))
      
      for(term.loop in colnames(shuffle.ts.z[[elec.loop]])){
        # term.loop = sp.model.terms[1]
        
        # Get z-score
        shuffle.ts.z[[elec.loop]][,term.loop] <-
          (current.shuffle.ts[[elec.loop]][,term.loop] - 
             shuffle.means[[elec.loop]][,term.loop]) /
          shuffle.sds[[elec.loop]][,term.loop]
        
        # Get sig windows at alpha = .05 for >= 100ms
        shuffle.zs.sig.windows[[elec.loop]][[term.loop]] <- list()
        for(significance.loop in names(significance.hps)){
          shuffle.zs.sig.windows[[elec.loop]][[term.loop]][[significance.loop]] <- 
            get.significant.windows(as.numeric(shuffle.ts.z[[elec.loop]][,term.loop] > significance.hps[[significance.loop]]['sig.threshold']),
                                    .sample.labels = rownames(shuffle.ts.z[[elec.loop]]),
                                    output.class = 'data.frame',
                                    include.duration = TRUE,
                                    .exclude.times.before.ms = time.convert(-median.rt.samples['sp'], "sample.labels", "times"),
                                    .exclude.times.after.ms = 250,
                                    .exclude.sig.durations.under.ms = significance.hps[[significance.loop]]['min.sig.window'])
        }; rm(significance.loop)
        
      }; rm(term.loop)
    }; rm(elec.loop)
    
    
    ### Get sig elecs
    shuffle.zs.sig.elecs <- list()
    output.n.sig.elecs.per.shuffle <- list()
    for(significance.loop in names(significance.hps)){
      
      shuffle.zs.sig.elecs[[significance.loop]] <- list()
      for(term.loop in names(shuffle.zs.sig.windows[[1]])){
        # term.loop = names(zs.sig.windows[[1]])[1]
        
        shuffle.zs.sig.elecs[[significance.loop]][[term.loop]] <- 
          sapply(shuffle.zs.sig.windows, function(x){
            nrow(x[[term.loop]][[significance.loop]]) > 0
          })
        shuffle.zs.sig.elecs[[significance.loop]][[term.loop]] <- 
          names(shuffle.zs.sig.elecs[[significance.loop]][[term.loop]][shuffle.zs.sig.elecs[[significance.loop]][[term.loop]]])
      }; rm(term.loop)
      
      # Get n.sig.elecs
      output.n.sig.elecs.per.shuffle[[significance.loop]] <- data.frame(t(sapply(shuffle.zs.sig.elecs[[significance.loop]], length)))
      
    }; rm(significance.loop)
    
    
    ### Adjust
    shuffle.zs.adjusted <- list()
    for(elec.loop in names(shuffle.ts.z)){
      # elec.loop = names(shuffle.ts.z)[1]
      
      # Storage
      shuffle.zs.adjusted[[elec.loop]] <- 
        data.frame(apply(shuffle.ts.z[[elec.loop]], 2, squish))
      
    }; rm(elec.loop)
    
    
    ###
    ### Save shuffle z's
    ### 
    
    save.these.shuffle.zs <- c(
      'shuffle.zs.adjusted',
      'shuffle.zs.sig.elecs'
    )
    save.shuffle.zs <- paste0(output.path, 'data/1c - stats - shuffled RSA/shuffle zs and sig elecs for each shuffle loop/')
    dir.create(save.shuffle.zs, showWarnings = FALSE, recursive = TRUE)
    save(list = save.these.shuffle.zs,
         file = paste0(save.shuffle.zs, 'zs and sig elecs - shuffle loop ',shuffle.loop,'.RData'))
    
    
    # Collect number of sig elecs per shuffle
    return(output.n.sig.elecs.per.shuffle)
    
  }#; rm(shuffle.loop)

## Close parallel backend
stopCluster(cl)
unregister_dopar()


### Plot n.sig.elecs vs. shuffle dist results
# Clean up
n.sig.elecs.in.shuffle.data <- list()
for(significance.loop in names(significance.hps)){
  n.sig.elecs.in.shuffle.data[[significance.loop]] <- 
    bind_rows(lapply(temp.n.sig.elecs.in.shuffle.data, function(x){x[[significance.loop]]}))
}; rm(significance.loop, temp.n.sig.elecs.in.shuffle.data)

# Plot directory
plot.n.sig.elecs.dir <- paste0(output.path, 'figures/1c - stats - shuffled RSA/number of sig elecs/')
dir.create(plot.n.sig.elecs.dir, showWarnings = FALSE, recursive = TRUE)

# Plot!
for(significance.loop in names(significance.hps)){
  # significance.loop = names(significance.hps)[4]
  
  # Get summary stats
  shuff.means <- apply(n.sig.elecs.in.shuffle.data[[significance.loop]], 2, mean)
  shuff.sds <- apply(n.sig.elecs.in.shuffle.data[[significance.loop]], 2, sd)
  
  pdf(paste0(plot.n.sig.elecs.dir, 'number of sig elecs - real data vs. shuffle distribution - ',significance.loop,'.pdf'),
       height = 4 * 3, width = 4)
  par(mfrow = c(3,1),
      oma = c(0,0,1,0))
  
  for(term.loop in colnames(n.sig.elecs[[significance.loop]])){
    # term.loop = colnames(n.sig.elecs[[significance.loop]])[1]
    current.z <- round((n.sig.elecs[[significance.loop]][,term.loop] - shuff.means[term.loop]) / shuff.sds[term.loop], 2)
    current.p <- round(1 - pnorm(current.z), 3)
    dist.from.real.to.mean <- abs(n.sig.elecs[[significance.loop]][,term.loop] - shuff.means[term.loop])
    hist(n.sig.elecs.in.shuffle.data[[significance.loop]][,term.loop], 
         bg = 'grey',
         main = paste0(term.loop,': z = ',current.z,'; p = ',current.p),
         bty = 'n',
         border = 'white',
         # breaks = 15,
         xlim = c(min(c(n.sig.elecs.in.shuffle.data[[significance.loop]][,term.loop], shuff.means[term.loop] - dist.from.real.to.mean)),
                  max(c(n.sig.elecs.in.shuffle.data[[significance.loop]][,term.loop], shuff.means[term.loop] + dist.from.real.to.mean))))
    abline(v = shuff.means[term.loop])
    segments(y0 = 0,
             x0 = shuff.means[term.loop] - (qnorm(.975) * shuff.sds[term.loop]),
             x1 = shuff.means[term.loop] + (qnorm(.975) * shuff.sds[term.loop]),
             col = 'red',
             lwd = 3)
    abline(v = n.sig.elecs[[significance.loop]][,term.loop],
           lwd = 2,
           col = 'royalblue')
  }; rm(term.loop)
  mtext(significance.loop,
        outer = TRUE,
        line = -1)
  dev.off()
}; rm(significance.loop)

### Save number of sig elecs per shuffle loop
save.n.shuffle.sig.elecs.path <- paste0(output.path, 'data/1c - stats - shuffled RSA/number of sig elecs/')
dir.create(save.n.shuffle.sig.elecs.path, showWarnings = FALSE, recursive = TRUE)
save(n.sig.elecs.in.shuffle.data,
     file = paste0(save.n.shuffle.sig.elecs.path, 'n sig elecs per shuffle loop.RData'))




# Finish!
message('Script completed successfully. ',Sys.time())
beep()






