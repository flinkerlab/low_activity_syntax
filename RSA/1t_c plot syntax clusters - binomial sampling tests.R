### Test regional makeup of clusters compared to each other and to original dataset distribution
### adam.milton.morgan@gmail.com


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
library('binom') # for binom.confint()
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

nmf.rank.loop = 3
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

###
### Plot individual clusters - high gamma and RSA terms, separately
###

### Load plot data (cluster.elec.stats, cluster.z.stats, cluster.elec.sig.windows, cluster.z.sig.windows)
load(paste0(output.path, 'data/1t - clustering just syntax/1t_b - values and stats for plots/',
            alpha.loop,'/',
            'rank=',nmf.rank.loop,
            '/syntax clusters - summary stats and plot values.RData'))

### Save peaks
save.data.dir <- paste0(output.path,
                        'data/1t - clustering just syntax/1t_c/RSI peaks by cluster/',
                        alpha.loop,'/',
                        'rank=',nmf.rank.loop,'/')
dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
cluster.rsi.peaks <- data.frame(bind_rows(lapply(cluster.z.stats, function(x){
  # x <- cluster.z.stats$NMF_1
  sapply(x, function(y){
    time.convert(y$sample.label[which.max(y$mean.wtd)], "sample.labels", "times")
  })}), .id = 'cluster'))
write.csv(cluster.rsi.peaks,
          paste0(save.data.dir, 'cluster RSI peaks - rank=',nmf.rank.loop,'.csv'),
          row.names = FALSE, quote = FALSE)

### Plot   
# Y-limits
ecog.ylims.max <- max(c(2, unlist(cluster.elec.maxes.weighted)))
rsa.ylims.max <- max(c(2, unlist(cluster.z.maxes.weighted)))
ecog.y.max.tick <- max(2, round(ecog.ylims.max))
rsa.y.max.tick <- max(2, round(rsa.ylims.max))

# Plot parameters
zoom <- 1.8

### Plots
for(theme.loop in c('black','white')){
  # theme.loop = 'white'
  
  term.colors <- c(setNames(colors[[theme.loop]]$rsa_term_line_plots$hex,
                            rownames(colors[[theme.loop]]$rsa_term_line_plots)),
                   c('diff.log.rt' = colors[[theme.loop]]$rainbow_bright['grey','hex']))
  save.singles.plot.path <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                   alpha.loop,'/',
                                   'individual models/ECoG and RSA time series/for publication/',
                                   'rank=',nmf.rank.loop,'/')
  dir.create(save.singles.plot.path, showWarnings = FALSE, recursive = TRUE)
  
  for(cluster.loop in sort(unique(cluster.assignments$cluster))){
    # cluster.loop = unique(cluster.assignments$cluster)[1]
    pdf(paste0(save.singles.plot.path, 'ecog - rank=', nmf.rank.loop,' - cluster=',cluster.loop,'.pdf'),
        width = 6, height = 3.5)
    par(oma = c(0,0,1,0))
    
    # ECoG in current band.loop
    task.order <- names(cluster.elec.stats[[cluster.loop]])
    plot.time.series(.y.values = lapply(cluster.elec.stats[[cluster.loop]][task.order], function(x){x$mean.wtd}),
                     .x.values = lapply(cluster.elec.stats[[cluster.loop]][task.order], function(x){x$sample.label}),
                     .error.bars = lapply(cluster.elec.stats[[cluster.loop]][task.order], function(x){x$se.wtd}),
                     # .sig.windows = cluster.elec.sig.windows[[cluster.loop]][task.order],
                     .x.limits = c(-max(median.rt.times) - 150, 500),
                     .y.limits = c(ifelse(band.loop == 'high_gamma', 0, -.25), ecog.ylims.max),
                     .y.ticks = c(0, ecog.y.max.tick),
                     .colors = colors[[theme.loop]]$task_line_plots[task.order,'hex'],
                     # .sig.color = colors[[theme.loop]]$task_line_plots[task.order,'hex'],
                     .theme = theme.loop,
                     .background = ifelse(theme.loop == 'white', rgb(1,1,1,0), rgb(0,0,0,1)),
                     # .horizontal.line.at = 0,
                     .y.label = '',
                     .x.ticks = c(-1000, -500, 0, 500),
                     .x.tick.labels = c('-1000', '', '0', '500'),
                     .x.label = '',
                     .zoom = zoom,
                     .margin = c(3,3,0,1))
    dev.off()
    
    
    # Z-scored t-values
    pdf(paste0(save.singles.plot.path, 'rsa - rank=', nmf.rank.loop,' - cluster=',cluster.loop,'.pdf'),
        width = 6, height = 3.5)
    par(oma = c(0,0,1,0))
    
    term.order <- names(cluster.z.stats[[cluster.loop]])
    plot.time.series(.y.values = lapply(cluster.z.stats[[cluster.loop]][term.order], function(x){x$mean.wtd}),
                     .x.values = lapply(cluster.z.stats[[cluster.loop]][term.order], function(x){x$sample.label}),
                     .error.bars = lapply(cluster.z.stats[[cluster.loop]][term.order], function(x){x$se.wtd}),
                     .sig.windows = cluster.z.sig.windows[[cluster.loop]][term.order],
                     .x.limits = c(-max(median.rt.times) - 150, 500),
                     .y.limits = c(0, rsa.ylims.max),
                     .y.ticks = c(0, rsa.y.max.tick),
                     .colors = term.colors[term.order],
                     .sig.color = term.colors[term.order],
                     # .horizontal.line.at = 0,
                     .theme = theme.loop,
                     .background = ifelse(theme.loop == 'white', rgb(1,1,1,0), rgb(0,0,0,1)),
                     .y.label = '',
                     .x.ticks = c(-1000, -500, 0, 500),
                     .x.tick.labels = c('-1000', '', '0', '500'),
                     .x.label = '',
                     .zoom = zoom,
                     .margin = c(3,3,0,1))
    
    dev.off()
    
  }; rm(cluster.loop)
  
}#; rm(theme.loop)




###
### Binomial tests
###

# Set up elec/region dataframe
cluster.regions.df <- cluster.assignments[,c('electrode','cluster')]
cluster.regions.df$region_clinical <- elec.info[cluster.regions.df$electrode,'region_clinical']
cluster.regions.df$roi <- gross.rois[cluster.regions.df$region_clinical,'region.category']
cluster.regions.df <- cluster.regions.df[which(!is.na(cluster.regions.df$roi)),] # get rid of uncategorized regions (should only be a few ~5ish)
cluster.regions.t <- with(cluster.regions.df, table(roi, cluster))

# Get useful counts etc.
cluster.regions.df <- as.data.frame.matrix(cluster.regions.t)[rois.to.plot,]
n.elecs.per.region.all <- table(all.roi.elecs.analyzed$gross.roi)
n.elecs.all <- sum(n.elecs.per.region.all)
proportion.elecs.per.region.all <- n.elecs.per.region.all / n.elecs.all
n.elecs.per.region.syntax <- table(all.roi.elecs.analyzed[all.roi.elecs.analyzed$syntax.sig == 1,]$gross.roi)
n.elecs.syntax <- sum(n.elecs.per.region.syntax)
rois <- rois.to.plot
n.rois <- length(rois)
n.elecs.per.cluster <- colSums(cluster.regions.df)
clusters <- sort(names(n.elecs.per.cluster))
n.clusters <- length(clusters)

### Binomial tests: is any region overrepresented relative to original data distribution?
cluster.roi.ps.all <- 
  cluster.roi.ps.syntax <-
  data.frame(matrix(nrow = nrow(cluster.regions.df),
                    ncol = ncol(cluster.regions.df),
                    dimnames = list(rownames(cluster.regions.df),
                                    colnames(cluster.regions.df))))
roi.ps.syntax <- data.frame(matrix(nrow = nrow(cluster.regions.df),
                                   ncol = 1,
                                   dimnames = list(rownames(cluster.regions.df),
                                                   'p')))
for(roi.loop in rois){
  # roi.loop = rois[1]
  
  # Distribution of syntax elecs (independent of cluster) compared to original distribution
  roi.ps.syntax[roi.loop, 'p'] <- 
    binom.test(x = n.elecs.per.region.syntax[roi.loop],
               n = n.elecs.syntax,
               p = n.elecs.per.region.all[roi.loop] / n.elecs.all,
               alternative = 'greater')$p.value
  
  for(cluster.loop in clusters){
    # cluster.loop = clusters[1]
    
    # Distribution of cluster elecs compared to original distribution
    cluster.roi.ps.all[roi.loop, cluster.loop] <-
      binom.test(x = cluster.regions.df[roi.loop, cluster.loop],
                 n = n.elecs.per.cluster[cluster.loop],
                 p = n.elecs.per.region.all[roi.loop] / n.elecs.all,
                 alternative = 'greater')$p.value
    # Distribution of cluster elecs compared to syntax distribution
    cluster.roi.ps.syntax[roi.loop, cluster.loop] <-
      binom.test(x = cluster.regions.df[roi.loop, cluster.loop],
                 n = n.elecs.per.cluster[cluster.loop],
                 p = n.elecs.per.region.syntax[roi.loop] / n.elecs.syntax,
                 alternative = 'greater')$p.value
    
  }; rm(cluster.loop)
}; rm(roi.loop)

# apply(cluster.roi.ps.all, c(1,2), function(x){as.numeric(x < .1)})
# apply(cluster.roi.ps.syntax, c(1,2), function(x){as.numeric(x < .1)})
# apply(roi.ps.syntax, c(1,2), function(x){as.numeric(x < .1)})

cluster.roi.ps.all.adj <- apply(cluster.roi.ps.all, 2, p.adjust, method = 'fdr')
cluster.roi.ps.syntax.adj <- apply(cluster.roi.ps.syntax, 2, p.adjust, method = 'fdr')
roi.ps.syntax.adj <- apply(roi.ps.syntax, 2, p.adjust, method = 'fdr')

# apply(cluster.roi.ps.all.adj, c(1,2), function(x){as.numeric(x < .05)})
# apply(cluster.roi.ps.syntax.adj, c(1,2), function(x){as.numeric(x < .1)})
# setNames(as.numeric(roi.ps.syntax.adj < .1), rownames(roi.ps.syntax.adj))


##
## Pie charts of region by cluster
##

# Hack the base r pie function to get rid of tick marks and color text the way you want it
my.pie <- function(x, 
                   labels = names(x),
                   .text.color = 'black',
                   .show.ticks = FALSE,
                   edges = 200,
                   radius = 0.8,
                   clockwise = TRUE,
                   init.angle = if (clockwise) 90 else 0,
                   density = NULL,
                   angle = 45,
                   col = NULL,
                   border = NULL,
                   .line.width = 1,
                   lty = NULL,
                   main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if(length(.text.color) < length(x)){
    .text.color <- rep(.text.color, length.out = length(x))
  }
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i], lwd = .line.width)
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      if(.show.ticks){lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)}
      text(1.1 * P$x,
           1.1 * P$y,
           labels[i],
           xpd = TRUE,
           col = .text.color[i],
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}

pie.plot.dir.count <- paste0(output.path, 'figures/1t - clustering just syntax/',
                             alpha.loop,'/',
                             'individual models/region-cluster pie charts/',
                             'rank=',nmf.rank.loop,
                             '/raw count/',
                             '/for publication/')
pie.plot.dir.proportion.all <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                      alpha.loop,'/',
                                      'individual models/region-cluster pie charts/',
                                      'rank=',nmf.rank.loop,
                                      '/proportion of all elecs/',
                                      '/for publication/')
pie.plot.dir.proportion.syntax <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                         alpha.loop,'/',
                                         'individual models/region-cluster pie charts/',
                                         'rank=',nmf.rank.loop,
                                         '/proportion of syntax elecs/',
                                         '/for publication/')
dir.create(pie.plot.dir.count, showWarnings = FALSE, recursive = TRUE)
dir.create(pie.plot.dir.proportion.all, showWarnings = FALSE, recursive = TRUE)
dir.create(pie.plot.dir.proportion.syntax, showWarnings = FALSE, recursive = TRUE)
for(theme.loop in c('white','black')){
  for(cluster.loop in clusters){
    # cluster.loop <- clusters[2]
    current.roi.labels <- rois.to.plot
    for(roi.loop in 1:length(current.roi.labels)){
      current.p <- round(cluster.roi.ps.all.adj[current.roi.labels[roi.loop], cluster.loop], 4)
      current.sig <- ""
      if(current.p < .1){current.sig <- "."}
      if(current.p < .05){current.sig <- "*"}
      if(current.p < .01){current.sig <- "**"}
      if(current.p < .001){current.sig <- "***"}
      current.roi.labels[roi.loop] <-
        paste0(current.roi.labels[roi.loop],
               ' p=',current.p,current.sig)
    }; rm(roi.loop)
    current.colors <- cubicl(n.rois)#[c(seq(1, n.rois, by = 2),seq(2, n.rois, by = 2))]
    
    ### Count pie plot
    pdf(paste0(pie.plot.dir.count,'cluster regions pie chart - ',theme.loop,' - rank=',nmf.rank.loop,' - ',cluster.loop,'.pdf'),
        height = 5, width = 5)
    par(oma = c(1,1,1,1),
        mar = c(0,0,0,0))
    my.pie(cluster.regions.df[rois.to.plot, cluster.loop],
           labels = "",
           .text.color = ifelse(theme.loop=="white","black","white"),
           border = theme.loop,
           col = current.colors,
           .line.width = 4)
    add.text.line.multiple.colors(current.roi.labels[1:3],
                                  text.colors = current.colors[1:3],
                                  .side = 3,
                                  .outer = FALSE)
    add.text.line.multiple.colors(current.roi.labels[4:6],
                                  text.colors = current.colors[4:6],
                                  .side = 1,
                                  .outer = FALSE)
    dev.off()
    
    ### Proportion all pie plot
    pdf(paste0(pie.plot.dir.proportion.all,'cluster regions pie chart - ',theme.loop,' - rank=',nmf.rank.loop,' - ',cluster.loop,'.pdf'),
        height = 5, width = 5)
    par(oma = c(1,1,1,1),
        mar = c(0,0,0,0))
    my.pie(cluster.regions.df[rois.to.plot, cluster.loop] / n.elecs.per.region.all[rois.to.plot],
           labels = "",
           .text.color = ifelse(theme.loop=="white","black","white"),
           border = theme.loop,
           col = current.colors,
           .line.width = 4)
    add.text.line.multiple.colors(current.roi.labels[1:3],
                                  text.colors = current.colors[1:3],
                                  .side = 3,
                                  .outer = FALSE)
    add.text.line.multiple.colors(current.roi.labels[4:6],
                                  text.colors = current.colors[4:6],
                                  .side = 1,
                                  .outer = FALSE)
    dev.off()
    
    ### Proportion syntax pie plot
    current.roi.labels <- rois.to.plot
    for(roi.loop in 1:length(current.roi.labels)){
      current.p <- round(cluster.roi.ps.syntax.adj[current.roi.labels[roi.loop], cluster.loop], 4)
      current.sig <- ""
      if(current.p < .1){current.sig <- "."}
      if(current.p < .05){current.sig <- "*"}
      if(current.p < .01){current.sig <- "**"}
      if(current.p < .001){current.sig <- "***"}
      current.roi.labels[roi.loop] <-
        paste0(current.roi.labels[roi.loop],
               ' p=',current.p,current.sig)
    }; rm(roi.loop)
    pdf(paste0(pie.plot.dir.proportion.syntax,'cluster regions pie chart - ',theme.loop,' - rank=',nmf.rank.loop,' - ',cluster.loop,'.pdf'),
        height = 5, width = 5)
    par(oma = c(1,1,1,1),
        mar = c(0,0,0,0))
    my.pie(cluster.regions.df[rois.to.plot, cluster.loop] / n.elecs.per.region.syntax[rois.to.plot],
           labels = "",
           .text.color = ifelse(theme.loop=="white","black","white"),
           border = theme.loop,
           col = current.colors,
           .line.width = 4)
    add.text.line.multiple.colors(current.roi.labels[1:3],
                                  text.colors = current.colors[1:3],
                                  .side = 3,
                                  .outer = FALSE)
    add.text.line.multiple.colors(current.roi.labels[4:6],
                                  text.colors = current.colors[4:6],
                                  .side = 1,
                                  .outer = FALSE)
    dev.off()
  }#; rm(cluster.loop)
}#; rm(theme.loop)


##
## Barplot
##

# Get proportions
cluster.roi.proportions <- cluster.regions.df
for(cluster.loop in colnames(cluster.roi.proportions)){
  cluster.roi.proportions[,cluster.loop] <-
    cluster.roi.proportions[,cluster.loop] / sum(cluster.roi.proportions[,cluster.loop])
}; rm(cluster.loop)

bar.plot.dir.clusters <- paste0(output.path, 'figures/1t - clustering just syntax/',
                                alpha.loop,'/',
                                'individual models/region-cluster bar charts/quick and dirty/by cluster/')
bar.plot.dir.rois <- paste0(output.path, 'figures/1t - clustering just syntax/',
                            alpha.loop,'/',
                            'individual models/region-cluster bar charts/quick and dirty/by roi/')
dir.create(bar.plot.dir.clusters, showWarnings = FALSE, recursive = TRUE)
dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)


# Add OG proportions
roi.proportions <- n.elecs.per.region.all / sum(n.elecs.per.region.all)
cluster.roi.proportions <-
  cbind(data.frame('all' = data.frame(roi.proportions[rownames(cluster.roi.proportions)])$Freq,
                   cluster.roi.proportions))

for(theme.loop in c('white','black')){
  
  # Plot by cluster
  roi.colors <- cubicl(n.rois)
  pdf(paste0(bar.plot.dir.clusters, 'cluster regions bar plot - ',theme.loop,' - rank=',nmf.rank.loop,'.pdf'),
      width = 1.8 * (1 + nmf.rank.loop), height = 4)
  par(oma = c(0,0,1,0), bg = rgb(1,1,1,1))
  barplot(data.matrix(cluster.roi.proportions),
          beside = TRUE,
          col = roi.colors,
          border = 'white')
  add.text.line.multiple.colors(rownames(cluster.roi.proportions),
                                text.colors = roi.colors,
                                .outer = TRUE)
  dev.off()
  
  
  # Plot by ROI
  # cluster.colors <- cubicl(n.clusters + 1)
  cluster.colors <- c(rgb(163, 0, 14, maxColorValue = 255), gray.colors(4)[1:3])
  pdf(paste0(bar.plot.dir.rois, 'cluster regions bar plot - ',theme.loop,' - rank=',nmf.rank.loop,'.pdf'),
      width = 1.8 * (1 + nmf.rank.loop), height = 4)
  par(oma = c(0,0,1,0), bg = rgb(1,1,1,1))
  barplot(t(data.matrix(cluster.roi.proportions)),
          beside = TRUE,
          col = cluster.colors,
          border = 'white')
  add.text.line.multiple.colors(colnames(cluster.roi.proportions),
                                text.colors = cluster.colors,
                                .outer = TRUE)
  dev.off()
  
  
  
  ### Plot by ROI cluster colored - FOR MANUSCRIPT
  
  ## Confidence intervals
  cis.lower <- cis.upper <- 
    data.frame(matrix(nrow = nrow(cluster.regions.df),
                      ncol = ncol(cluster.regions.df),
                      dimnames = list(rownames(cluster.regions.df),
                                      colnames(cluster.regions.df))))
  for(cluster.loop in colnames(cis.lower)){
    # cluster.loop = colnames(cis.lower)[3]
    for(roi.loop in rownames(cis.lower)){
      # roi.loop = rownames(cis.lower)[4]
      
      current.ci <- binom.confint(cluster.regions.df[roi.loop, cluster.loop],
                                  sum(cluster.regions.df[,cluster.loop]), 
                                  conf.level = .95, 
                                  methods = "wilson")
      cis.lower[roi.loop, cluster.loop] <- current.ci$lower
      cis.upper[roi.loop, cluster.loop] <- current.ci$upper
      
    }; rm(roi.loop)
  }; rm(cluster.loop)
  
  
  cluster.order <- c('NMF_2','NMF_3','NMF_1')
  cluster.colors <-
    c('NMF_2' = colors[[theme.loop]]$rainbow_lighter['orange','hex'],
      'NMF_3' = rgb(162,40,165, maxColorValue = 255), 
      'NMF_1' = rgb(40, 0, 80, maxColorValue = 255))
  
  
  
  ## Stats: one-vs-rest binomial tests
  roi.change.ps <- data.frame(matrix(nrow = length(rois.to.plot), 
                                     ncol = 3,
                                     dimnames = list(rois.to.plot, 
                                                     cluster.order)))
  
  for(cluster.loop in cluster.order){
    for(roi.loop in rois.to.plot){
      # More than expected based on other two clusters?
      current.test <-
        binom.test(x = cluster.regions.df[roi.loop, cluster.loop],
                   n = sum(cluster.regions.df[, cluster.loop]),
                   p = sum(cluster.regions.df[roi.loop, which(colnames(cluster.regions.df) != cluster.loop)]) /
                     sum(cluster.regions.df[, which(colnames(cluster.regions.df) != cluster.loop)]))
      roi.change.ps[roi.loop, cluster.loop] <- current.test$p.value
      
    }; rm(roi.loop)
  }; rm(cluster.loop) 
  
  roi.change.ps.adj <- apply(roi.change.ps, 2, p.adjust, method = 'fdr')
  
  # Text explaining what's significant
  sig.text <- "FDR-sig: "
  for(cluster.loop in clusters){
    # cluster.loop <- clusters[1]
    for(roi.loop in rois.to.plot){
      # roi.loop = rois.to.plot[1]
      current.p <- roi.change.ps.adj[roi.loop, cluster.loop]
      current.sig <- ""
      if(current.p < .1){current.sig <- "."}
      if(current.p < .05){current.sig <- "*"}
      if(current.p < .01){current.sig <- "**"}
      if(current.p < .001){current.sig <- "***"}
      if(current.sig != ""){
        sig.text <- paste0(sig.text, 
                           ' ',cluster.loop,'_',roi.loop,
                           ' p=',round(current.p, 4),current.sig,' & ')
      }
    }; rm(roi.loop)
  }; rm(cluster.loop)
}; rm(theme.loop)

text.size.med <- 1.6
zoom <- 1.8

for(theme.loop in c('black','white')){
  pdf(paste0(bar.plot.dir.rois, 'cluster regions bar plot - ',theme.loop,' - rank=',nmf.rank.loop,' - new.pdf'),
      width = 8.7, height = 6.95)
  par(oma = c(1,0,1,0),
      mar = c(4,4,0,0) * zoom,
      lwd = 2 * zoom,
      bg = rgb(1,1,1,0))
  barplot(t(data.matrix(cluster.roi.proportions[rois.to.plot, cluster.order])),
          ylim = c(0, max(cis.upper) * 1.02),
          beside = TRUE,
          col = cluster.colors[cluster.order],
          border = theme.loop,
          cex.names = text.size.med * zoom,
          yaxt = 'n',
          las = 2)
  # Error bars
  x <- .5
  for(roi.loop in rois.to.plot){
    # roi.loop = rois.to.plot[1]
    
    for(cluster.loop in cluster.order){
      # cluster.loop = cluster.order[1]
      
      x <- x + 1
      arrows(x0 = x,
             y0 = cis.lower[roi.loop, cluster.loop],
             y1 = cis.upper[roi.loop, cluster.loop],
             angle = 90,
             code = 3,
             length = 0, 
             col = rgb(.7,.7,.7))
    }; rm(cluster.loop)
    
    x <- x + 1
  }; rm(roi.loop)
  
  axis(side = 2,
       at = c(0, .5),
       # labels = .y.tick.labels,
       las = 0,
       tck = -.025 * zoom, # length of tick
       padj = -.45 * zoom, # distance between tick and label
       lwd = 1.5 * zoom,
       lwd.ticks = 1.5 * zoom,
       cex.axis = text.size.med * zoom,
       col = ifelse(theme.loop == 'white', 'black', 'white'),
       col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
  add.text.line.multiple.colors(text.segments = cluster.order,
                                text.colors = cluster.colors,
                                .side = 3,
                                .outer = TRUE)
  mtext(sig.text, 
        side = 1,
        outer = TRUE)
  dev.off()
  
}#; rm(theme.loop)


##
## Line charts of regions by cluster
##

for(theme.loop in rev(c('white','black'))){
  # theme.loop = 'white'
  pub.line.plot.dir <- paste0(output.path, 'figures/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/region-cluster line charts/counts/for publication/',
                              theme.loop,'/')
  dir.create(pub.line.plot.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - clusters.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = t(cluster.regions.df[rois.to.plot,]),
                  .p.values = t(cluster.roi.ps.all.adj[rois.to.plot,]),
                  .y.label = 'count',
                  show.legend = TRUE,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - regions.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = (cluster.regions.df[rois.to.plot,]),
                  .p.values = (cluster.roi.ps.all.adj[rois.to.plot,]),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  
  ### Proportion of elecs plots
  cluster.regions.proportions <- cluster.regions.df
  for(roi.loop in rois){
    cluster.regions.proportions[roi.loop,] <- cluster.regions.proportions[roi.loop,] / n.elecs.per.region.all[roi.loop]
  }; rm(roi.loop)
  for(cluster.loop in clusters){
    cluster.regions.proportions[,cluster.loop] <- cluster.regions.proportions[,cluster.loop] / n.elecs.per.cluster[cluster.loop]
  }; rm(cluster.loop)
  
  pub.line.plot.dir <- paste0(output.path, 'figures/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/region-cluster line charts/proportions/for publication/',
                              theme.loop,'/')
  dir.create(pub.line.plot.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - clusters.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = t(cluster.regions.proportions[rois.to.plot,]),
                  .p.values = t(cluster.roi.ps.all.adj[rois.to.plot,]),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - regions.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = (cluster.regions.proportions[rois.to.plot,]),
                  .p.values = (cluster.roi.ps.all.adj[rois.to.plot,]),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  
  ### PUBLICATION PLOTS
  
  
  ##
  ## Plot ratio of cluster proportion to OG proportion - line plots
  ##
  
  ### Ratio of cluster proportion to OG proportion - line plots
  cluster.regions.proportions <- cluster.regions.df
  for(cluster.loop in clusters){ # cluster.loop = clusters[1]
    cluster.regions.proportions[,cluster.loop] <- cluster.regions.proportions[,cluster.loop] / n.elecs.per.cluster[cluster.loop]
  }; rm(cluster.loop)
  for(roi.loop in rois){
    cluster.regions.proportions[roi.loop,] <- cluster.regions.proportions[roi.loop,] / proportion.elecs.per.region.all[roi.loop]
  }; rm(roi.loop)
  
  
  
  ### Convert CI counts to ratios for plotting
  roi.change.95ci.lower.ratio <- roi.change.95ci.lower.count
  roi.change.95ci.upper.ratio <- roi.change.95ci.upper.count
  for(cluster.loop in clusters){ # cluster.loop = clusters[1]
    roi.change.95ci.lower.ratio[,cluster.loop] <- roi.change.95ci.lower.ratio[,cluster.loop] / n.elecs.per.cluster[cluster.loop]
    roi.change.95ci.upper.ratio[,cluster.loop] <- roi.change.95ci.upper.ratio[,cluster.loop] / n.elecs.per.cluster[cluster.loop]
  }; rm(cluster.loop)
  for(roi.loop in rois){
    roi.change.95ci.lower.ratio[roi.loop,] <- roi.change.95ci.lower.ratio[roi.loop,] / proportion.elecs.per.region.all[roi.loop]
    roi.change.95ci.upper.ratio[roi.loop,] <- roi.change.95ci.upper.ratio[roi.loop,] / proportion.elecs.per.region.all[roi.loop]
  }; rm(roi.loop)
  
  
  ### Plot ratios (no 95% CIs because... it's really fucking complicated)
  pub.line.plot.dir <- paste0(output.path, 'figures/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/region-cluster line charts/ratio of OG and cluster proportions/for publication/',
                              theme.loop,'/')
  dir.create(pub.line.plot.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - clusters.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = t(cluster.regions.proportions[rois.to.plot,]),
                  .p.values = t(cluster.roi.ps.all.adj[rois.to.plot,]),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  # .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - regions.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = (cluster.regions.proportions[rois.to.plot,c('NMF_2','NMF_3','NMF_1')]),
                  # .p.values = (cluster.roi.ps.all.adj[rois.to.plot,]),
                  .colors = c(colors[[theme.loop]]$rainbow_lighter['orange','hex'],
                              colors[[theme.loop]]$rainbow_bright['purple','hex'],
                              rgb(0,0,0,1)),
                  .log.axis = "y",
                  .y.limits = c(.25,4),
                  # .y.limits = c(0,4),
                  .y.ticks = c(.25,.5,1,2,4),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  # .background = rgb(1,1,1,0),
                  # .error.bars.upper = (roi.change.95ci.upper.ratio[rois.to.plot,]),
                  # .error.bars.lower = (roi.change.95ci.lower.ratio[rois.to.plot,]),
                  .horizontal.line.at = 1,
                  .zoom = 1)
  dev.off()
  
  
  
  ### Plot proportions 
  pub.line.plot.dir <- paste0(output.path, 'figures/1t - clustering just syntax/',
                              alpha.loop,'/',
                              'individual models/region-cluster line charts/cluster proportions not scaled/for publication/',
                              theme.loop,'/')
  dir.create(pub.line.plot.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - clusters.pdf'),
      width = 8, height = 5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  plot.line.chart(.y.values = t(cluster.regions.proportions[rois.to.plot,]),
                  .p.values = t(cluster.roi.ps.all.adj[rois.to.plot,]),
                  show.legend = TRUE,
                  .theme = theme.loop,
                  # .background = rgb(1,1,1,0),
                  .zoom = 1)
  dev.off()
  
  
  
  
  #
  # Proportions within clusters
  #
  
  pdf(paste0(pub.line.plot.dir, 'shuffle - rank=', nmf.rank.loop,' - regions.pdf'),
      width = 10, height = 6.5)
  par(mfrow = c(1,1),
      oma = c(0,0,1,0))
  # Get proportions
  cluster.region.proportions <- list(
    'NMF_2' = cluster.regions.df[rois.to.plot, 'NMF_2'] / n.elecs.per.cluster['NMF_2'],
    'NMF_3' = cluster.regions.df[rois.to.plot, 'NMF_3'] / n.elecs.per.cluster['NMF_3'],
    'NMF_1' = cluster.regions.df[rois.to.plot, 'NMF_1'] / n.elecs.per.cluster['NMF_1']
  )
  plot.line.chart(.y.values = cluster.region.proportions,
                  .x.tick.labels = rois.to.plot,
                  # .p.values = (cluster.roi.ps.all.adj[rois.to.plot,]),
                  .colors = c(colors[[theme.loop]]$rainbow_lighter['orange','hex'],
                              colors[[theme.loop]]$rainbow_bright['purple','hex'],
                              rgb(0,0,0,1)),
                  # .y.limits = c(0,.6),
                  .jitter.x = TRUE,
                  .y.ticks = c(0, .5),
                  .x.left.pad = .25,
                  show.legend = TRUE,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .error.bars.upper = (roi.change.95ci.upper[rois.to.plot,]),
                  .error.bars.lower = (roi.change.95ci.lower[rois.to.plot,]),
                  .error.bar.lower.limit = 0,
                  .error.bar.upper.limit = 1,
                  .jitter.x.range = .24,
                  .zoom = 1.8)
  mtext(text = sig.text,
        side = 3,
        outer = TRUE,
        line = 0,
        col = 'black')
  dev.off()
  
  
  
} # theme.loop




# }; rm(nmf.rank.loop) # UNCOMMENT

# Clean up
# rm(list = save.these)
gc()

} # alpha.loop
# }; rm(band.loop)



message('Script completed successfully. ',Sys.time())
# Finish!









# } # include.hga.loop

