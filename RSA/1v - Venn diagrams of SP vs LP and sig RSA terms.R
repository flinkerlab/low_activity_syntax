### Plot Venn Diagrams (Euler plots) of overlap between (a) RSIs and (b) active/passive and sentence/list Wilcoxon tests
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
library('pals')
library('caret')
library('doParallel')
library('eulerr') # for Venn diagrams


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
    n.cores.to.use = 14
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 10 # big RDMs
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

### Colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

### RTs
load(paste0(path, 'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output/data/median RT samples.RData'))
median.rt.times <- setNames(time.convert(median.rt.samples, "samples", "times"), 
                            names(median.rt.samples))

### Loops
band.loop = c('high_gamma','beta')[1]

### Load sig elecs from RSA terms (RSIs)
load(paste0(path,"analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/",band.loop,"/data/1c - stats - shuffled RSA/real data stats/adjusted zs with sig windows and elecs.RData")) # Loads zs.adjusted, zs.sig.windows, zs.sig.elecs

#### Load sig elecs from active vs. passive Wilcoxon tests
load(paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/high_gamma/data/1r - syntax active-passive wilcoxon comparisons/sig windows and elecs/syntax sig windows and elecs.RData')) # Loads syntax.sig.elecs.df, syntax.wilcox.sig.elecs, syntax.wilcox.sig.windows

### Load sig elecs from sentence (SP) vs. list (LP) Wilcoxon tests
load(paste0(path,'analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/data/warped data/high_gamma/patients - combined/from stimulus to 400ms post speech onset/',
            'sp vs. lp comparisons - all times to and from median rt of sp.RData')) # Loads sp.greater.lp.sig.elecs, sp.greater.lp.sig.windows, sp.greater.than.lp.ps
# 'sp vs. lp comparisons - all times to and from median rt of lp.RData'))

# Metadata
term.labels <- c('diff.voice' = 'syntax',
                 'diff.1st.word' = 'word',
                 'diff.event.semantics' = 'semantics')
sig.hps <- intersect(names(sp.greater.lp.sig.elecs), intersect(names(zs.sig.elecs), names(syntax.wilcox.sig.elecs)))


##
## Re-define sig elecs to just those sig during planning period (from stimulus to speech)
##

### zs.sig.elecs (RSIs)
# Update sig windows
for(elec.loop in names(zs.sig.windows)){
  for(term.loop in names(zs.sig.windows[[elec.loop]])){
    for(sig.hp.loop in sig.hps){
      temp <- zs.sig.windows[[elec.loop]][[term.loop]][[sig.hp.loop]]
      # If any sig windows...
      if(nrow(temp) > 0){
        # ...then for each sig window
        for(i in 1:nrow(temp)){
          # Set the earliest and latest times to stim onset and speech onset
          temp$start.time[i] <- max(temp$start.time[i], -median.rt.times['sp'])
          temp$end.time[i] <- min(temp$end.time[i], 0)
        }; rm(i)
        # Re-calculate durations
        temp$duration <- temp$end.time - temp$start.time
        # Filter
        temp <- temp[temp$duration >= as.numeric(gsub("alpha=...+_min=","",gsub("ms","",sig.hps))),]
        # Store
        zs.sig.windows[[elec.loop]][[term.loop]][[sig.hp.loop]] <- temp
        rm(temp)
      } # if(nrow(temp) > 0){
    }; rm(sig.hp.loop)
  }; rm(term.loop)
}; rm(elec.loop)

# Redefine zs.sig.elecs
zs.sig.elecs <- list()
for(sig.hp.loop in sig.hps){
  zs.sig.elecs[[sig.hp.loop]] <- list()
  for(term.loop in names(zs.adjusted[[1]])){
    zs.sig.elecs[[sig.hp.loop]][[term.loop]] <- 
      names(which(sapply(zs.sig.windows, function(x){
        nrow(x[[term.loop]][[sig.hp.loop]]) > 0
      })))
  }; rm(term.loop)
}; rm(sig.hp.loop)

### Do the same for syntax.wilcox.sig.windows and sp.greater.lp.sig.windows (structured the same)
for(sig.hp.loop in sig.hps){
  # sig.hp.loop = sig.hps[1]
  syntax.wilcox.sig.windows[[sig.hp.loop]] <- 
    lapply(syntax.wilcox.sig.windows[[sig.hp.loop]], function(temp){
      if(nrow(temp) > 0){
        # ...then for each sig window
        for(i in 1:nrow(temp)){
          # Set the earliest and latest times to stim onset and speech onset
          temp$start.time[i] <- max(temp$start.time[i], -median.rt.times['sp'])
          temp$end.time[i] <- min(temp$end.time[i], 0)
        }; rm(i)
        # Re-calculate durations
        temp$duration <- temp$end.time - temp$start.time
        # Filter
        temp <- temp[temp$duration >= as.numeric(gsub("alpha=...+_min=","",gsub("ms","",sig.hps))),]
      } # if(nrow(temp) > 0){
      return(temp)
    })
  sp.greater.lp.sig.windows[[sig.hp.loop]] <- 
    lapply(sp.greater.lp.sig.windows[[sig.hp.loop]], function(temp){
      if(nrow(temp) > 0){
        # ...then for each sig window
        for(i in 1:nrow(temp)){
          # Set the earliest and latest times to stim onset and speech onset
          temp$start.time[i] <- max(temp$start.time[i], -median.rt.times['sp'])
          temp$end.time[i] <- min(temp$end.time[i], 0)
        }; rm(i)
        # Re-calculate durations
        temp$duration <- temp$end.time - temp$start.time
        # Filter
        temp <- temp[temp$duration >= as.numeric(gsub("alpha=...+_min=","",gsub("ms","",sig.hps))),]
      } # if(nrow(temp) > 0){
      return(temp)
    })
}; rm(sig.hp.loop)

# Redefine sig.elecs
syntax.wilcox.sig.elecs <- 
  lapply(syntax.wilcox.sig.windows, function(x){
    names(which(sapply(x, function(y){
      nrow(y) > 0
    })))
  })
sp.greater.lp.sig.elecs <-
  lapply(sp.greater.lp.sig.windows, function(x){
    names(which(sapply(x, function(y){
      nrow(y) > 0
    })))
  })


# Only look at elecs from patients who produced passives -- otherwise unfairly biases SP>LP comparison
patients.who.produced.passives <- unique(substr(zs.sig.elecs$`alpha=.05_min=100ms`$diff.voice, 1, 5))
use.these.elecs <- use.these.elecs[substr(use.these.elecs, 1, 5) %in% patients.who.produced.passives]

# Set up data template
elec.df.template <- data.frame('all' = rep(TRUE, length(use.these.elecs)),
                               row.names = use.these.elecs)

# for(sig.hp.loop in sig.hps){ ## UNCOMMENT
sig.hp.loop = sig.hps[1]

###
### Venn Diagrams (Euler Plots)
###


### Set up comparison of RSA results and ECoG (SP>LP) results
ecog.rsa.df <- elec.df.template
ecog.rsa.df$ecog <- rownames(ecog.rsa.df) %in% sp.greater.lp.sig.elecs[[sig.hp.loop]]
for(term.loop in names(term.labels)){
  # term.loop = names(term.labels)[1]
  ecog.rsa.df[,term.loop] <- rownames(ecog.rsa.df) %in% zs.sig.elecs[[sig.hp.loop]][[term.loop]]
} # term.loop

### Set up comparison of ECOG (SP>LP) and Wilocxon (Active != Passive) results
ecog.cox.df <- elec.df.template
ecog.cox.df$ecog <- rownames(ecog.cox.df) %in% sp.greater.lp.sig.elecs[[sig.hp.loop]]
ecog.cox.df$diff.voice <- rownames(ecog.cox.df) %in% syntax.wilcox.sig.elecs[[sig.hp.loop]]


### Save plots
circle.transparency <- .65
for(theme.loop in c('white','black')){ ## UNCOMMENT
  # theme.loop = 'white'
  
  all.elecs.color <- colors[[theme.loop]]$rainbow_bright['cyan','hex']
  ecog.sp.color <- colors[[theme.loop]]$task_line_plots['sp','hex']
  
  ### ECoG, RSA syntax, RSA semantics
  current.plot.label <- 'ecog and RSA - syntax, semantics'
  current.terms <- c('diff.voice', 'diff.event.semantics')
  current.colors <- c(ecog.sp.color, colors[[theme.loop]]$rsa_term_line_plots[c(current.terms),'hex'])
  save.venn.dir <- paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/figures/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',theme.loop,'/',current.plot.label,'/')
  dir.create(save.venn.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - with labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   quantities = TRUE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - no labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   labels = c('', rep('', length(current.terms))),
                   quantities = FALSE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  
  ### ECoG, RSA syntax, RSA semantics, all
  current.plot.label <- 'all, ecog, and RSA - syntax, semantics'
  current.terms <- c('diff.voice', 'diff.event.semantics')
  current.colors <- c(all.elecs.color, ecog.sp.color, colors[[theme.loop]]$rsa_term_line_plots[c(current.terms),'hex'])
  save.venn.dir <- paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/figures/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',theme.loop,'/',current.plot.label,'/')
  dir.create(save.venn.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - with labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('all', 'ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   quantities = TRUE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - no labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('all','ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   labels = c('', '', rep('', length(current.terms))),
                   quantities = FALSE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  
  ### RSA syntax, RSA semantics, RSA word
  current.plot.label <- 'RSA - syntax, semantics, word'
  current.terms <- c('diff.voice', 'diff.event.semantics', 'diff.1st.word')
  current.colors <- c(colors[[theme.loop]]$rsa_term_line_plots[c(current.terms),'hex'])
  save.venn.dir <- paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/figures/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',theme.loop,'/',current.plot.label,'/')
  dir.create(save.venn.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - with labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c(current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   quantities = TRUE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - no labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c(current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   labels = c('', '', rep('', length(current.terms))),
                   quantities = FALSE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  
  ### ECoG and RSA syntax
  current.plot.label <- 'ecog and RSA - syntax'
  current.terms <- c('diff.voice')
  current.colors <- c(ecog.sp.color, colors[[theme.loop]]$rsa_term_line_plots[c(current.terms),'hex'])
  save.venn.dir <- paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/figures/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',theme.loop,'/',current.plot.label,'/')
  dir.create(save.venn.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - with labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   quantities = TRUE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - no labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.rsa.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   labels = c('', rep('', length(current.terms))),
                   quantities = FALSE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  
  ### ECoG and Wilcoxon active/pasive comaprison syntax
  current.plot.label <- 'ecog and Wilcoxon syntax'
  current.terms <- c('diff.voice')
  current.colors <- c(ecog.sp.color, colors[[theme.loop]]$rsa_term_line_plots[c(current.terms),'hex'])
  save.venn.dir <- paste0(path,'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',band.loop,'/figures/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',theme.loop,'/',current.plot.label,'/')
  dir.create(save.venn.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - with labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.cox.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   quantities = TRUE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  pdf(paste0(save.venn.dir, current.plot.label,' - ',sig.hp.loop,' - no labels.pdf'),
      width = 5,height = 5)
  fit <- euler(ecog.cox.df[,c('ecog', current.terms)])
  fit.plot <- plot(fit,
                   fill = current.colors,
                   edges = FALSE,
                   alpha = circle.transparency,
                   labels = c('', rep('', length(current.terms))),
                   quantities = FALSE)
  par(pty = 's',
      bg = rgb(1,1,1,0))
  print(fit.plot)
  dev.off()
  
  
}#; rm(theme.loop)


###
### Save data for brain plots
###

# Make numeric
for(col.loop in 1:ncol(ecog.cox.df)){
  if(class(ecog.cox.df[,col.loop]) == "logical"){ecog.cox.df[,col.loop] <- as.numeric(ecog.cox.df[,col.loop])}
}

# Add localizations etc
ecog.cox.df$elec <- rownames(ecog.cox.df)
ecog.cox.df <- cbind(ecog.cox.df,
                     elec.info[ecog.cox.df$elec, c('region_clinical', 'MNI_x','MNI_y','MNI_z')])
colnames(ecog.cox.df) <- gsub(".", "_", colnames(ecog.cox.df), fixed = TRUE)
save.data.dir <- paste0(path,
                        'analysis/R/event and syntactic encoding/elecs/rsa/rsa encoding models - no CV - clustering just SP ts - stats from shuffling/output - step 1/',
                        band.loop,
                        '/data/1v - Venn diagrams of SP vs LP wilcoxon tests and RSA terms/',
                        sig.hp.loop,'/')
dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
write.csv(ecog.cox.df,
          paste0(save.data.dir, 'sentence vs list and active vs passive sig elecs.csv'),
          row.names = FALSE, quote = FALSE)


# }; rm(sig.hp.loop)



























message('Script successfully completed. ', Sys.time())









