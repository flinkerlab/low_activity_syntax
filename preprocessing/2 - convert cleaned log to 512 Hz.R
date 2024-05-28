### Chicken Syntax
### Convert parsed psychopy log to 512 Hz log
### May 2020
### adam.milton.morgan@gamil.com

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
  path = '/Users/am4611/Dropbox/Research/ChickenSyntax/data/'} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'} # Ubuntu

### Read data
log <- read.csv(paste0(path,patient,'/log/parsed/',patient,'_cleaned_log.csv'))
# Make new 512 Hz dataframe
log.long <- data.frame('time'=time <- seq(0,max(log$time),by=1/512))

# Make new columns
log.long$diode.on <- 0
log.long$image.on <- 0
log.long$text.on <- 0

log.long$task.type <- NA
log.long$image.character <- NA
log.long$image.character.orientation <- NA
log.long$image.action <- NA
log.long$image.active.subject.position <- NA
log.long$image.id <- NA
log.long$image.arrow <- NA
log.long$text <- NA
log.long$active.or.passive <- NA
log.long$quick.image.duration <- NA

log.long$phase <- NA
log.long$task <- NA
log.long$onset.trial <- 0
log.long$onset.phase <- 0
log.long$onset.task.name.character <- 0
log.long$onset.task.read.question <- 0
log.long$onset.task.describe.scene <- 0
log.long$onset.task.read.sentence <- 0
log.long$onset.task.see.arrow <- 0
log.long$onset.task.list.characters <- 0
log.long$onset.task.read.list <- 0

# Useful variable
nrow.log.long <- nrow(log.long)

# Initialize big loop
# Strategy: Going row by row thru the event ("short") log,
# Every time there's a change, change everything in log long from that time onward
# So every new event overwrites the previous event from the new event start time onward
(start.time <- Sys.time())
for(row in 1:nrow(log)){ # row <- 4608
  # Find relevant rows in log.long
  log.long.row <- which(abs(log.long$time - log$time[row])==
                          min(abs(log.long$time - log$time[row])))
  log.long$phase[log.long.row:nrow.log.long] <- as.character(log$phase[row])
  if(log$event_diode[row]==1){
    if(as.character(log$command[row])=="on"){
      log.long$diode.on[log.long.row:nrow.log.long] <- 1
    }
    if(as.character(log$command[row])=="off"){
      log.long$diode.on[log.long.row:nrow.log.long] <- 0
    }
  }
  if(!is.na(log$task_type[row])){
    log.long$task.type[log.long.row:nrow.log.long] <- as.character(log$task_type[row])
  }
  if(log$event_text[row]==1 & (! is.na(log$command[row]))){
    if(as.character(log$command[row])=="on"){
      log.long$text.on[log.long.row:nrow.log.long] <- 1
    }
    if(as.character(log$command[row])=="off"){
      log.long$text.on[log.long.row:nrow.log.long] <- 0
    }
  }
  if(log$event_text[row]==1 & ! is.na(log$command[row])){
    if(as.character(log$command[row])=="on"){
      log.long$text[log.long.row:nrow.log.long] <- as.character(log$text_id[row])
    }
    if(as.character(log$command[row])=="off"){
      log.long$text[log.long.row:nrow.log.long] <- NA
    }
  }
  ## Add whether the question prompt is active or passive voice
  if(!is.na(log$active_or_passive[row])){
    log.long$active.or.passive[log.long.row:nrow.log.long] <- as.character(log$active_or_passive[row])
  }
  # End it
  if(! log$task_type[row] %in% c('produce_sentence','read_sentence','read_question',NA)){
    log.long$active.or.passive[log.long.row:nrow.log.long] <- NA
  }
  ## Image stuff
  if(log$event_image[row]==1 & (! is.na(log$command[row])) & log$event_id[row]!="diode_on"){
    if(as.character(log$command[row])=="on"){
      log.long$image.on[log.long.row:nrow.log.long] <- 1
      log.long$image.character[log.long.row:nrow.log.long] <- as.character(log$image_character_id[row])
      log.long$image.character.orientation[log.long.row:nrow.log.long] <- as.character(log$image_character_orientation[row])
      log.long$image.action[log.long.row:nrow.log.long] <- as.character(log$image_action[row])
      log.long$image.active.subject.position[log.long.row:nrow.log.long] <- as.character(log$image_active_subject_position[row])
      log.long$image.id[log.long.row:nrow.log.long] <- as.character(log$image_id[row])
      if(! is.na(log$image_arrow[row])){
        log.long$image.arrow[log.long.row:nrow.log.long] <- as.character(log$image_arrow[row])  
      }
    }
    if(as.character(log$command[row])=="off"){
      log.long$image.on[log.long.row:nrow.log.long] <- 0
      log.long$image.character[log.long.row:nrow.log.long] <- NA
      log.long$image.character.orientation[log.long.row:nrow.log.long] <- NA
      log.long$image.action[log.long.row:nrow.log.long] <- NA
      log.long$image.active.subject.position[log.long.row:nrow.log.long] <- NA
      # log.long$image.id[log.long.row:nrow.log.long] <- NA
      # log.long$image.arrow[log.long.row:nrow.log.long] <- NA
    }
  }
  ## Keep track of stimulus information throughout the trial, even after stimulus image is off screen:
  # Arrows
  if(! log$task_type[row] %in% c('see_arrow','produce_list',NA)){
    log.long$image.arrow[log.long.row:nrow.log.long] <- NA
  }
  # Scenes
  if(! log$task_type[row] %in% c('produce_sentence','read_sentence','produce_list',NA)){
    log.long$image.id[log.long.row:nrow.log.long] <- NA
  }
  if(log$new_phase[row]==1){log.long$onset.phase[log.long.row] <- 1}
  if(log$new_trial[row]==1){
    log.long$onset.trial[log.long.row] <- 1
    # Add image duration if quick naming phase
    if(log$phase[row]=="character_naming_quick"){
      log.long$quick.image.duration[log.long.row:nrow.log.long] <- log$image_duration[row]
    }
  }
  if((! is.na(log$task_type[row])) & (! is.na(log$command[row]))){
    if(as.character(log$command[row])=="on"){
      log.long$task[log.long.row:nrow.log.long] <- as.character(log$task_type[row])
      if(as.character(log$task_type[row]) == 'name_character'){
        log.long$onset.task.name.character[log.long.row] <- 1
      }else{
        if(log$phase[row] %in% c('sentence_production','sentence_production_again')){
          if(log$task_type[row]=='read_question'){
            log.long$onset.task.read.question[log.long.row] <- 1}
          if(log$task_type[row]=='produce_sentence'){
            log.long$onset.task.describe.scene[log.long.row] <- 1}
          if(log$task_type[row]=='read_sentence'){
            log.long$onset.task.read.sentence[log.long.row] <- 1}
        }
        if(log$phase[row]=='list_production'){
          if(log$task_type[row]=='see_arrow'){
            log.long$onset.task.see.arrow[log.long.row] <- 1}
          if(log$task_type[row]=='produce_list'){
            log.long$onset.task.list.characters[log.long.row] <- 1}
          if(log$task_type[row]=='read_list'){
            log.long$onset.task.read.list[log.long.row] <- 1}
          }
        }
      }
    if(as.character(log$command[row])=="off"){
      log.long$task[log.long.row:nrow.log.long] <- NA
    }
  }

  if((row%%500)==0){message("Time to complete row ",row," of ",nrow(log)," (",round(100*row/nrow(log),2),"%): ",round((Sys.time()-start.time),2))}
}
message("Total big loop duration: ",round(Sys.time()-start.time,2)) # 7 mins

# Replace periods in column names with underscores
colnames(log.long) <- gsub(".","_",colnames(log.long),fixed=TRUE)

# Picture is only on the screen for 300ms in "sentence_production_again" phase, so task is listed as NA after this. fix that
log.long$task[which(log.long$phase == 'sentence_production_again' & is.na(log.long$task))] <- 'produce_sentence'

# Set image duration to NA for all phases where images are on screen for unlimited time
log.long$quick_image_duration[log.long$phase != "character_naming_quick"] <- NA

# Factorize
log.long$task_type <- as.factor(log.long$task_type)
log.long$image_character <- as.factor(log.long$image_character)
log.long$image_character_orientation <- as.factor(log.long$image_character_orientation)
log.long$image_action <- as.factor(log.long$image_action)
log.long$image_active_subject_position <- as.factor(log.long$image_active_subject_position)
log.long$image_arrow <- as.factor(log.long$image_arrow)
log.long$text <- as.factor(log.long$text)
log.long$active_or_passive <- as.factor(log.long$active_or_passive)
log.long$phase <- as.factor(log.long$phase)
log.long$task <- as.factor(log.long$task)

# Write
write.csv(log.long,paste0(path,patient,'/log/parsed/',patient,'_log_512_no_words.csv'),row.names=FALSE,quote=FALSE)

