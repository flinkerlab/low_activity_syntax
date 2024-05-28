### Figure out time offsets between ECoG, audio, and PsychoPy log reconstructions
### May 2020
### adam.milton.morgan@gmail.com

rm(list=ls())
cat("\014")

#####################
### Manually edit ###
#####################
patient <- 'Patient010'

#############################################
### Auto - more manual edits needed below ###
#############################################

### Set path
if(Sys.info()['sysname'] == 'Darwin'){
  path = '/Users/am4611/Dropbox/Research/ChickenSyntax/data/'} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'} # Ubuntu

### Lemmas
library('tuneR')
library('signal')

###
### Figure out offset between audio and ECoG
###

### Load data
ecog_audio <- read.csv(paste0(path,patient,'/data/',patient,'_ECoG_audio_channel.txt'), header=FALSE, col.names="ecog_audio")$ecog_audio
audio_512 <- read.csv(paste0(path,patient,'/audio/',patient,'_audio_512.txt'), header=FALSE, col.names="audio")$audio
word_log <- read.csv(paste0(path,patient,'/transcription/',patient,'_word_log.csv'))

# Visualize
plot.start <- 30000
plot.end <- 60000
plot(ecog_audio[plot.start:plot.end]/max(abs(ecog_audio[plot.start:plot.end]))+.5, type='l', ylim=c(-1.5,1.5), col='coral')
lines(audio_512[plot.start:plot.end]/max(abs(audio_512[plot.start:plot.end]))/4 - .5, type='l', col='aquamarine')
text(x = 10000, y = c(.1,-.1), labels = c('ECoG (DC)','Audio (Recorder)'), col = c('coral','aquamarine'), cex = 1.2)

#####################
### Manually edit ###
#####################

# Looks like audio begins and ends before ECoG, offset by 80k-ish samples
# Clip off or add to beginning of audio:
samples.removed.from.beginning.of.audio <- 0
samples.added.to.beginning.of.audio <- 32400

############
### Auto ###
############

audio_512_adjusted <- c(rep(0, times=samples.added.to.beginning.of.audio),
                        audio_512)[(samples.removed.from.beginning.of.audio+1):length(audio_512)]
lines(audio_512_adjusted[plot.start:plot.end]/max(abs(audio_512[plot.start:plot.end]))/4 + .5, type='l', col='purple')


# Cross correlate audio and ECoG_audio
lag.max <- 2000
# # Old approach: cross-correlate the whole thing. BUT...
# xc <- ccf(ecog_audio, audio_512_adjusted, type="correlation", lag.max=lag.max)
# # ... since there's drift, we want to find out the offset at the beginning and end and then try to correct
xc <- ccf(ecog_audio[plot.start:plot.end], audio_512_adjusted[plot.start:plot.end], type="correlation", lag.max=lag.max)
(relative.audio.offset <- which(abs(xc$acf)==max(abs(xc$acf))) - lag.max) # +63 for Patient010

# Do the offset
if(relative.audio.offset < 0){
  samples.removed.from.beginning.of.audio <- samples.removed.from.beginning.of.audio - relative.audio.offset # minus a negative
}
if(relative.audio.offset > 0){
  samples.added.to.beginning.of.audio <- samples.added.to.beginning.of.audio + relative.audio.offset
}

# Perform the correct adjustment
audio_512_adjusted <- c(rep(0, times=samples.added.to.beginning.of.audio),
                        audio_512)[(samples.removed.from.beginning.of.audio+1):length(audio_512)]

########
######## BEGIN: Updated to include drift correciton, paste into other Preprocessing 4 scripts
########

#####################
### Manually edit ###
#####################

plot.end <- 100000
# Check visually
ecog.short <- ecog_audio[plot.start:plot.end]
audio.short <- audio_512_adjusted[plot.start:plot.end]
plot((ecog.short - mean(ecog.short))/max(abs(ecog.short)),
     col='coral',
     type='l',lty=1)
lines((audio.short - mean(audio.short))/max(abs(audio.short))/8,
      col='purple', #'aquamarine',#rgb(127,252,212,alpha=100,maxColorValue=400),
      type='l',lty=1)

# Check again: zoom in a little
sub.start <- 35000
sub.end <- 55000
ecog.short <- ecog_audio[sub.start:sub.end]
audio.short <- audio_512_adjusted[sub.start:sub.end]
plot((ecog.short - mean(ecog.short))/max(abs(ecog.short)),
     col='coral',
     type='l',lty=1)
lines((audio.short - mean(audio.short))/max(abs(audio.short))/6,
      col='purple',
      type='l',lty=1)

## Check the drift by looking at the end
sub.end <- min(length(ecog_audio), length(audio_512_adjusted))
sub.end <- sub.end - 10000
sub.start <- sub.end - 20000

ecog.short <- ecog_audio[sub.start:sub.end]
audio.short <- audio_512_adjusted[sub.start:sub.end]
plot((ecog.short - mean(ecog.short))/max(abs(ecog.short)),
     col='coral',
     type='l',lty=1)
lines((audio.short - mean(audio.short))/max(abs(audio.short))/6,
      col='purple',
      type='l',lty=1)

# Zoom in a little
sub.start.zoom <- sub.start
sub.end.zoom <- sub.start + 10000
ecog.short.zoom <- ecog_audio[sub.start.zoom:sub.end.zoom]
audio.short.zoom <- audio_512_adjusted[sub.start.zoom:sub.end.zoom]
plot((ecog.short.zoom - mean(ecog.short.zoom))/max(abs(ecog.short.zoom)),
     col='coral',
     type='l',lty=1)
lines((audio.short.zoom - mean(audio.short.zoom))/max(abs(audio.short.zoom))/6 - .1,
      col='purple',
      type='l',lty=1)

# How offset?
drift.xc <- ccf(ecog_audio[sub.start:sub.end], audio_512_adjusted[sub.start:sub.end], type="correlation", lag.max=lag.max)
(drift.relative.audio.offset <- which(abs(drift.xc$acf)==max(abs(drift.xc$acf))) - lag.max) # +51 for Patient010

## Visualize correcting for the offset at the end (due to drift now):
# Perform the adjustment on just the subset of data to plot
audio.short.adjusted <- 
  c(rep(0, times=ifelse(drift.relative.audio.offset > 0, drift.relative.audio.offset, 0)),
                          audio.short)[(1 + ifelse(drift.relative.audio.offset < 0, abs(drift.relative.audio.offset), 0)):length(audio.short)]

# Visualize
plot((ecog.short - mean(ecog.short))/max(abs(ecog.short)),
     col='coral',
     type='l',lty=1)
lines((audio.short.adjusted - mean(audio.short.adjusted))/max(abs(audio.short.adjusted))/8,
      col='purple',
      type='l',lty=1)
# Looks aligned! 

# So the drift rate is:
# NY 857: 43 samples / (1670000 - 10000) = 2.590361e-05 samples per sample, with an intercept at 10000 (plot.start)
(drift.rate <- drift.relative.audio.offset / (sub.start - plot.start))

# If it's a linear drift, we should be able to calculate how far off the recording is at any given moment:
drift.calculator <- function(sample = sub.start, round = TRUE){
  output <- (sample - plot.start) * drift.rate
  if(round){output <- round(output, 0)}
  return(output)
}

# Let's check to see if that works (it may not: the drift could be non-linear)
test.drift.calc.start <- min(length(ecog_audio), length(audio_512_adjusted)) - 20000
test.drift.calc.end <- test.drift.calc.start + 20000
test.drift.audio.short <- audio_512_adjusted[test.drift.calc.start:test.drift.calc.end]
test.drift.ecog.short <- ecog_audio[test.drift.calc.start:test.drift.calc.end]
(test.drift.calc.drift <- drift.calculator(sample = test.drift.calc.start))
test.drift.calc.audio.short.adjusted <- 
  c(rep(0, times=ifelse(test.drift.calc.drift > 0, test.drift.calc.drift, 0)),
    test.drift.audio.short)[(1 + ifelse(test.drift.calc.drift < 0, abs(test.drift.calc.drift), 0)):length(audio.short)]
plot((test.drift.ecog.short - mean(test.drift.ecog.short))/max(abs(test.drift.ecog.short)),
     col='coral',
     type='l',lty=1)
# OG (with drift) in purple
lines((test.drift.audio.short - mean(test.drift.audio.short))/max(abs(test.drift.audio.short))/6,
      col='purple',
      type='l',lty=1)
# Corrected in aquamarine
lines((test.drift.calc.audio.short.adjusted - mean(test.drift.calc.audio.short.adjusted))/max(abs(test.drift.calc.audio.short.adjusted))/6,
      col='green4',
      type='l',lty=1)
# Zoom:
plot((test.drift.ecog.short[6000:8500] - mean(test.drift.ecog.short))/max(abs(test.drift.ecog.short)),
     col='coral',
     type='l',lty=1)
lines((test.drift.audio.short[6000:8500] - mean(test.drift.audio.short))/max(abs(test.drift.audio.short))/6,
      col='purple',
      type='l',lty=1)
lines((test.drift.calc.audio.short.adjusted[6000:8500] - mean(test.drift.calc.audio.short.adjusted))/max(abs(test.drift.calc.audio.short.adjusted))/6,
      col='green4',
      type='l',lty=1)
abline(v = c(445, 525, 1780), lwd = .5, lty = 1)

############
### Auto ###
############

## Get number of seconds to add to or remove from beginning of audio
(seconds.added.to.beginning.of.audio <- round((samples.added.to.beginning.of.audio - samples.removed.from.beginning.of.audio) / 512, 2))  # should be negative if removing time from beginning of audio

# Remove rows with empty start/end times
word_log <- droplevels(subset(word_log, !((begin_s == "") | is.na(begin_s) | (end_s=="") | is.na(end_s))))

word_log$begin_s <- word_log$begin_s + seconds.added.to.beginning.of.audio
word_log$end_s <- word_log$end_s + seconds.added.to.beginning.of.audio


## Adjust for drift
# Convert to samples
word_log$begin_samples <- word_log$begin_s * 512
word_log$end_samples <- word_log$end_s * 512

# Correct
word_log$begin_samples <- word_log$begin_samples + drift.calculator(sample = word_log$begin_samples, round = FALSE)
word_log$end_samples <- word_log$end_samples + drift.calculator(sample = word_log$end_samples, round = FALSE)

# Save as times again
word_log$begin_s <- word_log$begin_samples / 512
word_log$end_s <- word_log$end_samples / 512

# Clean up
word_log$begin_samples <- NULL
word_log$end_samples <- NULL


# Write
write.csv(word_log,
          paste0(path,
                 patient,
                 '/transcription/',
                 patient,
                 '_word_log_aligned_w_ECoG.csv'),
          row.names=FALSE,
          quote=FALSE)

########
######## END: Updated to include drift correciton, paste into other Preprocessing 4 scripts
########



###
### Figure out offset between log and ECoG
###

rm(list=ls()[which(!ls() %in% c('path','patient'))])
log <- read.csv(paste0(path,patient,'/log/parsed/',patient,'_log_512_no_words.csv'))
diode <- read.table(paste0(path,patient,'/data/',patient,'_photodiode_channel.txt'), header=FALSE)$V1
beep <- read.table(paste0(path,patient,'/data/',patient,'_trigger_channel.txt'), header=FALSE)$V1

#####################
### Manually edit ###
#####################

## Look at just a snippet of the beginning of the task
short.begin.row <- 4e4#50000
short.end.row <- 150000

############
### Auto ###
############

# Plot
diode.short <- diode[short.begin.row:short.end.row]/max(abs(diode[short.begin.row:short.end.row]))
beep.short <- beep[short.begin.row:short.end.row]/max(abs(beep[short.begin.row:short.end.row]))
log.short <- log[(short.begin.row):(short.end.row),]
plot(1-diode.short+1.2, type='l', lty=1,ylim=c(-1,3))
lines(log.short$diode_on+.5, type='l',lty=1,col='blue')
lines(beep.short, type='l',lty=1,col='red')


#####################
### Manually edit ###
#####################

samples.removed.from.beginning.of.log <- 7.5e4

short.begin.row <- 2e4
short.end.row <- 8e4

############
### Auto ###
############

diode.short <- diode[short.begin.row:short.end.row]/max(abs(diode[short.begin.row:short.end.row]))
beep.short <- beep[short.begin.row:short.end.row]/max(abs(beep[short.begin.row:short.end.row]))
log.short <- log[(short.begin.row+samples.removed.from.beginning.of.log):(short.end.row+samples.removed.from.beginning.of.log),]
plot(1-diode.short+1.2, type='l', lty=1,ylim=c(-1,3))
lines(log.short$diode_on+.5, type='l',lty=1,col='blue')
lines(beep.short, type='l',lty=1,col='red')

# Cross-correlate
lag.max <- 1000
xc <- ccf(diode.short, log.short$diode_on, type="correlation", lag.max=lag.max)
(relative.diode.offset <- which(abs(xc$acf)==max(abs(xc$acf))) - lag.max) 

# Do offset
samples.removed.from.beginning.of.log <- samples.removed.from.beginning.of.log - relative.diode.offset
log.short <- log[(short.begin.row+samples.removed.from.beginning.of.log):(short.end.row+samples.removed.from.beginning.of.log),]

# Check again
lag.max <- 1000
xc <- ccf(diode.short, log.short$diode_on, type="correlation", lag.max=lag.max)
which(abs(xc$acf)==max(abs(xc$acf))) - lag.max

# Check plot again
short.begin.row <- 500000
short.end.row <- 650000
diode.short <- diode[short.begin.row:short.end.row]/max(abs(diode[short.begin.row:short.end.row]))
beep.short <- beep[short.begin.row:short.end.row]/max(abs(beep[short.begin.row:short.end.row]))
log.short <- log[(short.begin.row+samples.removed.from.beginning.of.log):(short.end.row+samples.removed.from.beginning.of.log),]
plot(1-diode.short+1.2, type='l', lty=1,ylim=c(-1,3))
lines(log.short$diode_on+.5, type='l',lty=1,col='blue')
lines(beep.short, type='l',lty=1,col='red')

# How many samples to remove?
samples.removed.from.beginning.of.log # for NY829: 18035, NY834: 2202

# Write aligned log
log.aligned <- log[samples.removed.from.beginning.of.log:nrow(log),]
log.aligned$sample <- 1:nrow(log.aligned)
(subtract.time <- log.aligned$time[1])
log.aligned$time <- log.aligned$time - subtract.time

write.csv(log.aligned,
          paste0(path,
                 patient,
                 '/log/parsed/',
                 patient,
                 '_log_512_aligned_w_ECoG.csv'),
          row.names=FALSE,
          quote=FALSE)

