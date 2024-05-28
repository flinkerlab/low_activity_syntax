### Define MFG axis, model electrodes' peak syntax RSI time as function of projected position along MFG axis
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
# Adjust color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Create sigmoids
source(paste0(path,'/analysis/R/functions/sigmoid.R'))
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))


## RTs
load(paste0(path,'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
median.rt.times <- setNames(time.convert(median.rt.samples, "samples", "times"), names(median.rt.samples))

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
               (elec.info$bad_elec == 0))
rownames(elec.info) <- elec.info$patient_elec

use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec  

## Save directory
output.path <- paste0(path, 'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/')


### Load electrode importances
elec.zs.path <- paste0(output.path,'data/1c - stats - shuffled RSA/real data stats/')
elec.zs.file <- 'adjusted zs with sig windows and elecs.RData'
load(paste0(elec.zs.path, elec.zs.file))


# Significance hyperparameters to loop through
alpha.hps <- c('alpha=.05_min=100ms',
               'alpha=.01_min=50ms')

# for(alpha.loop in alpha.hps){
alpha.loop = alpha.hps[1]

# NMF ranks to try
nmf.ranks.to.try <- 3

roi.clinical.subregions <- list(
  'mfg' = c('rostralmiddlefrontal','caudalmiddlefrontal')
)
rois = names(roi.clinical.subregions)
# for(roi.loop in rois){
roi.loop = rois[1]

## Read term times
term.times <- read.csv(paste0(output.path,'data/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/csvs for brain plots/rank=',
                              nmf.ranks.to.try,
                              '/term_times.csv'))

# Clean up
term.times <- term.times[,c('electrode','MNI_y','MNI_z','region_clinical','peak_time_diff_voice','cluster')]
term.times <- term.times[term.times$region_clinical %in% roi.clinical.subregions[[roi.loop]],]

### Get a projected MFG axis
# Perform PCA on the x and y coordinates
pca_result <- prcomp(term.times[, c("MNI_y", "MNI_z")], scale. = FALSE)

# Get the unit vector of the first principal component
principal_axis <- pca_result$rotation[, 1]

# Project electrode coordinates onto the principal axis
projections <- as.matrix(term.times[, c("MNI_y", "MNI_z")]) %*% principal_axis

# Add the projections as a new column in the dataframe
term.times$projection <- projections

# Get slope and angle of projected axis
slope <- principal_axis[2] / principal_axis[1]
intercept <- mean(term.times$MNI_z) - slope * mean(term.times$MNI_y)
angle <- atan(slope) * 180 / pi # in degrees

# Plot the principal axis to visualize this new dimension
plot(y = term.times$MNI_z,
     x = term.times$MNI_y,
     ylim = rev(range(term.times$MNI_z)))

# Add the principal axis line to the plot
abline(a = intercept, b = slope, col = "red", lwd = 2)

### Limit to just data in the planning period
term.times <- term.times[term.times$peak_time_diff_voice <= 0,]
term.times <- term.times[term.times$peak_time_diff_voice > -median.rt.times['sp'],]


## Save dir
save.dir <- paste0(output.path,
                   'figures/1t - clustering just syntax/1t_d - by roi/scatterplots - MNI by time/PC1 of MNI_y and MNI_z/')
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

# Plot hyperparameters
text.size.big <- 1.8
text.size.med <- 1.6
text.size.small <- 1.4
zoom <- 1.8

## Term colored plot
model <- with(term.times, lm(peak_time_diff_voice ~ projection))
current.p <- summary(model)$coefficients[2,4]

for(theme.loop in c('black','white')){
  # theme.loop = 'white'
  
  # Load colors
  colors.lightest <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',theme.loop,'/rainbow_lightest.csv'), row.names = 1)
  colors.lighter <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',theme.loop,'/rainbow_lighter.csv'), row.names = 1)
  colors.bright <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',theme.loop,'/rainbow_bright.csv'), row.names = 1)
  colors.darker <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',theme.loop,'/rainbow_darker.csv'), row.names = 1)
  colors.darkest <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',theme.loop,'/rainbow_darkest.csv'), row.names = 1)
  
  # Add cluster-specific colors
  cluster.colors <- 
    c('NMF_2' = colors.lighter['orange','hex'],
      'NMF_3' = rgb(162,40,165, maxColorValue = 255), 
      'NMF_1' = rgb(40, 0, 80, maxColorValue = 255))
  # c('NMF_2' = colors.lighter['orange','hex'], # early
  #   'NMF_3' = colors.bright['purple','hex'], # middle
  #   'NMF_1' = rgb(0,0,0,1))
  # setNames(rev(magma(8)[c(1,4,7)]), c('NMF_2','NMF_3','NMF_1'))
  term.times$cluster.color <- cluster.colors[term.times$cluster]
  
  # Add subregion (rostral/caudal) specific colors
  subroi.colors <- 
    c('rostralmiddlefrontal' = colors.lighter['orange','hex'], # early
      'caudalmiddlefrontal' = colors.lighter['purple','hex']) # late
  term.times$subroi.color <- subroi.colors[term.times$region_clinical]
  
  for(color.code.loop in c('cluster.color','subroi.color')){
    
    ### Plot
    
    # Begin save
    pdf(paste0(save.dir, 'MNI (projection of Y-Z plane) by peak syntax time - ',
               color.code.loop,' - ',
               theme.loop,
               ' - p=',round(current.p, 3),
               ' - angle of projection=',round(angle, 2),'.pdf'),
        width = 5.7, height = 5.7)
    
    # Plot set up
    time.range <- c(-median.rt.times['sp'], 0)
    MNI.range <- range(term.times$projection)
    MNI.range <- ceiling(abs(MNI.range)/5) * 5 * sign(MNI.range)
    par(bg = ifelse(theme.loop == 'black', rgb(0,0,0,0), rgb(1,1,1,0)),
        # mar = c(6,6,3,1) * zoom,
        mar = c(3,3,0,1) * zoom,
        oma = c(0,0,1,0))
    
    # Plot
    plot(x = term.times$projection,
         y = term.times$peak_time_diff_voice,
         xlab = '',
         xlim = rev(MNI.range + c(0, 2)),
         ylim = time.range + c(-50, 0),
         ylab = '',
         pch = 20,
         cex = 3 * zoom,
         col = adjust.transparency(term.times[,color.code.loop], 1),
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n')
    abline(model,
           lwd = 4 * zoom,
           # col = rgb(.5,.5,.5))
           col = ifelse(theme.loop == 'black', 'white', 'black'))
    time.at <- c(-median.rt.times['sp'], 0)
    time.at.labels <- as.character(time.at)
    MNI.at <- MNI.range
    MNI.at.labels <- as.character(MNI.at)
    axis(side = 1,
         at = MNI.at,
         labels = MNI.at.labels,
         las = 1,
         tck = -.025 * zoom, # length of tick
         padj = .6 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.med * zoom,
         col = ifelse(theme.loop == 'black', 'white', 'black'),
         col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
    axis(side = 2,
         at = time.at,
         labels = time.at.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.med * zoom,
         col = ifelse(theme.loop == 'black', 'white', 'black'),
         col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
    
    dev.off()
    
    
    
    
    
  }#; rm(color.code.loop)           
}#; rm(theme.loop)
# }#; rm(roi.loop)          
#         } # alpha.loop
# # }; rm(band.loop)



# Finish!
message('Script completed successfully. ',Sys.time())









