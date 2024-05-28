### Chicken syntax
### Convert .log file to something legible to then turn into a 512 sampling rate
### adam.milton.morgan@gmail.com
### May 2020

rm(list=ls())

#####################
### Manually edit ###
#####################
patient <- 'Patient010'

####################################
### No manual edits needed below ###
####################################

### Set path
if(Sys.info()['sysname'] == 'Darwin'){
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/data/'} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'} # Ubuntu

### Packages
library(splitstackshape)
library(dplyr)
library(stringr)
library(tidyr)

### Read in data
log.filename <- list.files(paste0(path,patient,'/log/log_file/'))
log <- read.table(paste0(path,patient,'/log/log_file/',log.filename),
                  sep='\t',
                  col.names=c('time','del','output'))
log <- droplevels(subset(log, del=='EXP '))
log$del <- NULL

# Separate output column into various columns
log <- separate(log, 'output', into=c('name','command_og','del','rest'), sep=" ", extra="merge")
# Get rid of beginning rows, etc.
log <- droplevels(subset(log, command_og %in% c('trial','text','image','autoDraw') &
                           name != 'Created'))
log$del <- NULL
log$name <- as.factor(log$name)
log$command_og <- as.factor(log$command_og)

# Get "quick" time out of rest column
log$image_duration <- NA
log$image_duration[grepl("(quick_duration, 300)", log$rest)] <- 300
log$image_duration[grepl("(quick_duration, 200)", log$rest)] <- 200
log$rest[which(log$command_og=='trial')] <- NA

## Phase
log$phase <- NA
log$new_phase <- 0
phases <- c('intro',
            'character_naming',
            'sentence_production',
            'list_production',
            "character_naming_again",
            "character_naming_ecologically_valid",
            "character_naming_quick",
            "sentence_production_again")
log$phase[1] <- phases[1]
j <- 1
new.phase.markers <- c("pre_train_text:",
                       "practice_scene_text:",
                       "pre_order_text:",
                       "pre_name_text:",
                       "pre_eco_name_text:",
                       "pre_naming_last_text:",
                       "pre_description_last_text:")
for(i in 2:nrow(log)){
  if(log$name[i] %in% new.phase.markers & log$rest[i]=="True"){
    j <- j+1
    log$new_phase[i] <- 1
  }
  log$phase[i] <- phases[j]
}; rm(i,j, phases, new.phase.markers)
log$phase <- as.factor(log$phase)

# Set image_duration to NA for all phases except character_naming_quick
log$image_duration[log$phase != "character_naming_quick"] <- NA

## Trial
log$trial <- NA
log$new_trial <- 0
log$trial[1] <- 0
for(i in 2:nrow(log)){
  if(log$command_og[i]=="trial"){
    log$trial[i] <- log$trial[i-1]+1
    log$new_trial[i] <- 1
  }else{
    log$trial[i] <- log$trial[i-1]
    }
}; rm(i)

## Event type (diode, image, beep, text, new_trial, new_phase)
log$event_diode <- 0
log$event_image <- 0
log$event_text <- 0
# log$event_sound <- 0
for(i in 1:nrow(log)){
  temp.name <- as.character(log$name[i])
  if(temp.name != gsub("diode","",temp.name)){
    log$event_diode[i] <- 1
  }
  if(temp.name != gsub("order_prime_arrow","",gsub("image","",temp.name))){
    log$event_image[i] <- 1
  }
  if((temp.name != gsub("text","",temp.name) | temp.name == "scene_A:") &
     (! temp.name %in% c("text:","null_text:","order_prime_text:"))){
    log$event_text[i] <- 1
  }
}; rm(i,temp.name)

## Event command (on, off, NA)
log$command <- NA
log$command[with(log, which(command_og=="autoDraw" & rest=="True" & (event_diode==1 | event_image==1 | event_text==1)))] <- "on"
log$command[with(log, which(command_og=="autoDraw" & rest=="False" & (event_diode==1 | event_image==1 | event_text==1)))] <- "off"
log$command <- as.factor(log$command)

## Event content (image_id, text_id)
# Image ID column
log$image_id <- NA
for(i in 1:nrow(log)){
  if(log$event_image[i]==1){
    if(as.character(log$command_og[i])=="image"){
      temp.pic <- temp.pic <- gsub("scenes/","",gsub("characters/","",gsub(".png","",gsub("images/","",as.character(log$rest[i])))))
    }
    if(log$command_og[i]=="autoDraw" & log$command[i]=="on"){
      log$image_id[i] <- temp.pic
    }
    if(log$command_og[i]=="autoDraw" & log$command[i]=="off"){
      log$image_id[i] <- temp.pic
      temp.pic <- NA
    }
  }
}; rm(i,temp.pic)
log$image_id <- as.factor(log$image_id)

# Text ID column
log$text_id <- NA
temp.text <- NA
for(i in 1:nrow(log)){
  if(log$event_text[i]==1){
    if(log$command_og[i]=="text"){
      temp.text <- as.character(log$rest[i])
    }
    if(log$command_og[i]=="autoDraw" & log$command[i]=="on"){
      log$text_id[i] <- temp.text
    }
    if(log$command_og[i]=="autoDraw" & log$command[i]=="off"){
      log$text_id[i] <- temp.text
      temp.text <- NA
    }
  }
}; rm(i,temp.text)
log$text_id <- as.factor(log$text_id)

# Clean up
log <- droplevels(subset(log, new_phase==1 | new_trial==1 | event_diode==1 | event_image==1 | event_text==1))
log <- droplevels(subset(log, command_og %in% c('autoDraw','trial')))
log$rest <- NULL
log$command_og <- NULL
row.names(log) <- NULL

# Add task_type column (name_char, read_q, read_sentence, read_list, see_arrow)
log$task_type <- NA
log$task_type[which(log$name %in% c('char1_image:','char2_image:','train_image:','name_image:','name_quick_image:','name_ecological_image:','name_last_image:'))] <- "name_character"
log$task_type[which(log$name %in% c('scene_image:','scene_last_image:'))] <- 'produce_sentence'
log$task_type[which(log$name %in% c('order_image:'))] <- 'produce_list'
log$task_type[which(log$name %in% c('scene_prime_text:'))] <- "read_question"
log$task_type[which(log$name %in% c('scene_prime_verb_text:'))] <- "read_verb"
log$task_type[which(log$name %in% c('scene_A:'))] <- "read_sentence"
log$task_type[which(log$name %in% c('order_A_text:'))] <- "read_list"
log$task_type[which(log$name %in% c('order_prime_arrow:'))] <- "see_arrow"
log$task_type <- as.factor(log$task_type)

# Add active or passive column
log$active_or_passive <- NA
for(i in 1:nrow(log)){
  if(as.character(log$task_type[i]) %in% c('read_question','produce_sentence','read_sentence')){
    if(as.character(log$task_type[i]) == 'read_question'){
      if(as.character(log$text_id[i]) != gsub(" is being ","",gsub(" is getting ","",as.character(log$text_id[i])))){
        temp.act.pass <- 'passive'
      }else{temp.act.pass <- 'active'}
    }
    if(log$command[i]=="on"){
      log$active_or_passive[i] <- temp.act.pass
    }
    if(log$command[i]=="off"){
      log$active_or_passive[i] <- temp.act.pass
      temp.pic <- NA
    }
  }
}; rm(i, temp.act.pass)

# Make the last sentence phase all NAs for act/pass since the prompt is just a verb
log$active_or_passive[which(log$phase=='sentence_production_again')] <- NA

# Add event_id column
log$event_id <- NA
log$event_id[which(log$new_phase==1)] <- 'new_phase'
log$event_id[which(log$new_trial==1)] <- 'new_trial'
log$event_id[which(log$event_diode==1 & log$command=='on')] <- 'diode_on'
log$event_id[which(log$event_diode==1 & log$command=='off')] <- 'diode_off'
log$event_id[which(log$event_image==1 & log$command=='on')] <- 'image_on'
log$event_id[which(log$event_image==1 & log$command=='off')] <- 'image_off'
log$event_id[which(log$event_text==1 & log$command=='on')] <- 'text_on'
log$event_id[which(log$event_text==1 & log$command=='off')] <- 'text_off'
log$event_id <- as.factor(log$event_id)

# Add image character ID columns for name_caracter task
# Arrow images
log$image_arrow <- NA
log$image_arrow[which(log$task_type=='see_arrow')] <- as.character(log$image_id[which(log$task_type=='see_arrow')])
log$image_arrow <- as.factor(log$image_arrow)
log$image_id[which(log$task_type=='see_arrow')] <- NA

# Split image_id column up
log <- cSplit(log, 'image_id', sep="-", type.convert=TRUE, drop=FALSE)
log$image_id <- as.factor(log$image_id)
# Add image_character_id column
log$image_character_id <- NA
log$image_character_id[which(log$task_type=="name_character")] <- as.character(log$image_id_1[which(log$task_type=="name_character")])
log$image_character_id <- as.factor(log$image_character_id)

# Add image_character_orientation column
log$image_character_orientation <- NA
log$image_character_orientation[which(log$task_type=="name_character")] <- as.character(log$image_id_2[which(log$task_type=="name_character")])
log$image_character_orientation <- as.numeric(log$image_character_orientation)

# Add image_action column
log$image_action <- NA
log$image_action[which(log$task_type %in% c("produce_sentence","produce_list"))] <- as.character(log$image_id_1[which(log$task_type %in% c("produce_sentence","produce_list"))])
log$image_action <- as.factor(log$image_action)

# Add image_active_subject_position column
log$image_active_subject_position <- NA
log$image_active_subject_position[which(log$task_type %in% c("produce_sentence","produce_list"))] <- as.character(log$image_id_3[which(log$task_type %in% c("produce_sentence","produce_list"))])
log$image_active_subject_position <- as.factor(log$image_active_subject_position)

# Add image_condition column
log$image_condition <- NA
log$image_condition[which(log$task_type %in% c("produce_sentence","produce_list"))] <- as.character(log$image_id_2[which(log$task_type %in% c("produce_sentence","produce_list"))])
log$image_condition <- as.factor(log$image_condition)

# Clean up
log$image_id_1 <- NULL
log$image_id_2 <- NULL
log$image_id_3 <- NULL

# Write
if(! dir.exists(paste0(path,patient,'/log/parsed/'))){dir.create(paste0(path,patient,'/log/parsed/'))}
write.csv(log, file=paste0(path,patient,'/log/parsed/',patient,'_cleaned_log.csv'), row.names=FALSE, quote=FALSE)

