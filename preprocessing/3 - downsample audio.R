### Chicken Syntax
### Downsample audio and save various triggers
### May 2020
### adam.milton.morgan@gamil.com

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
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/data/'} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'} # Ubuntu

### Packages, functions, data, etc.
library('edfReader')
library('tuneR')
library('signal')
library('beepr')
center <- function(x){
  x <- (x - mean(x)) / sd(x)
}

### Get data
# Read in data
header = readEdfHeader(paste0(path,patient,'/data/',patient,'_ChickenSyntax_512.EDF'))
data = readEdfSignals(header)
# Read in audio
audio = readWave(paste0(path,patient,'/audio/',patient,'_ChickenSyntax_uni_mono.wav'))
samp.rate = audio@samp.rate
beep()

####################
### Manual check ###
####################

# Confirm that the audio is mono from the unidirectional mic visually
plot(extractWave(audio, from=800000, to=1200000))

############
### Auto ###
############

# Downsample (use left channel even if mono)
audio.512 = resample(audio@left, p=128, q=11025)

# Write downsampled audio
write(audio.512,
      paste0(path,patient,'/audio/',patient,'_audio_512.txt'),
      ncolumns=1)


#########################################################
### Manually figure out which channel is which signal ###
#########################################################

# Plot audio and diode triggers and audio recroding
plot(center(data$DC5$signal[700000:710000]), type='l', col='lightblue') # noise
lines(center(data$DC2$signal[700000:710000]), col='red') # audio trigger
lines(center(data$DC3$signal[700000:710000]) / 5, col='darkgreen') # audio recording
lines(center(data$DC4$signal[700000:710000]), col='blue') # noise
lines(center(data$DC1$signal[700000:710000]), col='black', lwd=3) # diode

# 24/7: noisy audio, picking up audio in whole room # no longer collected for patients after 7/2021
# (photo)diode: square tooth with values of 0 or 1, delimits trials
# audio trigger: values -1, 0, 1; denotes various events; more common than photo diode
# recorder out: another noisey audio recording, can use this or 24/7 to align with high quality recording

trigger_data <- data$DC2$signal # audio trigger channel
diode_data <- data$DC1$signal # photo diode channel
ecog_audio_data <- data$DC3$signal # 24/7 or recorder out, whichever's cleanest

#######################
### Auto - the rest ###
#######################

# Write channel data - different channels for different participants!!
write(trigger_data, 
      paste0(path,patient,'/data/',patient,'_trigger_channel.txt'),
      ncolumns=1)
write(diode_data, 
      paste0(path,patient,'/data/',patient,'_photodiode_channel.txt'),
      ncolumns=1)
write(ecog_audio_data, 
      paste0(path,patient,'/data/',patient,'_ECoG_audio_channel.txt'),
      ncolumns=1)

