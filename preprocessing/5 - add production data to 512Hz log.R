### Chicken Syntax
### Add production data to 512 Hz log
### May 2020
### adam.milton.morgan@gmail.com

rm(list=ls())

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

# Read data
log <- read.csv(paste0(path,patient,'/transcription/',patient,'_word_log_aligned_w_ECoG.csv'))
log.long <- read.csv(paste0(path,patient,'/log/parsed/',patient,'_log_512_aligned_w_ECoG.csv'))

# Clean up log
# Get rid of extra rows (may need these later at some point)
log$comment <- NULL

# Remove any words uttered in less than .02s
log$duration <- with(log, end_s - begin_s)
if(length(which(log$duration<=0))>0){log <- log[-which(log$duration<=0),]} # drop any rows with 0 or negative word durations

#####################
### Manually edit ###
#####################

log$word <- as.character(log$word)
log$word <- tolower(log$word)
sort(unique(log$word))
log$word[which(log$word=="don't")] <- 'dont'
log$word[which(log$word=="scar")] <- 'scare'
log$word[which(log$word=="frankenstien")] <- 'frankenstein'
log$word[which(log$word=="ok")] <- 'okay'
log$word[which(log$word=="don\xd5t")] <- 'dont'
log$word[which(log$word=="dont")] <- 'dont'
log$word[which(log$word=="that\xd5s")] <- 'thats'
# sort(summary(as.factor(log$word)))
table(log$word)
# sort(summary(as.factor(droplevels(subset(log, pos=='n'))$word)))


##########################################
### Auto - no more manual edits  below ###
##########################################

# Get rid of rows without words we care about (6 nouns) & "exclude" rows
log <- droplevels(subset(log, pos != "" & 
                           (exclude_from_syntax_analyses==0 | 
                              is.na(exclude_from_syntax_analyses)) & 
                           ! is.na(duration)))
# Clean up
log$pos <- as.character(log$pos)
log$dp_or_np <- as.character(log$dp_or_np)
log$role <- as.character(log$role)

# Act/pass congruent column
if(ncol(with(log[log$incong_produced != "",], table(incong_produced, behavior))) > 1){
  message("Warning! Incong column off:")
  print(with(log[log$incong_produced != "",], table(incong_produced, behavior)))
  message("Go back and fix in .xslx -> CSV -> 2 -> 4 -> 5")
}
log$incong_produced[(log$behavior == "sentence") & (log$incong_produced == "")] <- "cong"
log$incong_produced[log$behavior != "sentence"] <- NA

# Make new columns
log.long$word <- NA
log.long$onset.chicken <- 0
log.long$onset.dog <- 0
log.long$onset.dracula <- 0
log.long$onset.frankenstein <- 0
log.long$onset.ninja <- 0
log.long$onset.nurse <- 0
log.long$pos <- NA
log.long$dp_or_np <- NA
log.long$syntactic_role <- NA
log.long$onset_aux <- 0
log.long$onset_verb <- 0
log.long$onset_determiner <- 0
log.long$onset_by <- 0
log.long$onset_passive_be <- 0
log.long$noun1 <- NA
log.long$noun2 <- NA
log.long$verb <- NA
log.long$verb_voice <- NA
### BEGIN ADDED 10/21/2022
log.long$substitution_error <- NA
log.long$correct_noun <- NA
log.long$wrong_noun <- NA
### END ADDED 10/21/2022
log.long$voice_congruent <- NA


# Store common variables
nrow.log.long <- nrow(log.long)

# Initialize big loop
# Strategy: Going row by row thru the event ("short") log,
# Every time there's a change, change everything in log long from that time onward
# So every new event overwrites the previous event from the new event start time onward
(start.time <- Sys.time())
for(row in 1:nrow(log)){ 
  # row <- 31 # for troubleshooting
  # Find relevant rows in log.long
  log.long.start.row <- which(abs(log.long$time - log$begin_s[row])==
                                min(abs(log.long$time - log$begin_s[row])))[1]
  log.long.end.row <- which(abs(log.long$time - log$end_s[row])==
                              min(abs(log.long$time - log$end_s[row])))[1]
  # Add in new variables
  w <- as.character(log$word[row])
  log.long$word[log.long.start.row:log.long.end.row] <- w
  if(w=='chicken'){log.long$onset.chicken[log.long.start.row] <- 1}
  if(w=='dog'){log.long$onset.dog[log.long.start.row] <- 1}
  if(w=='dracula'){log.long$onset.dracula[log.long.start.row] <- 1}
  if(w=='frankenstein'){log.long$onset.frankenstein[log.long.start.row] <- 1}
  if(w=='ninja'){log.long$onset.ninja[log.long.start.row] <- 1}
  if(w=='nurse'){log.long$onset.nurse[log.long.start.row] <- 1}
  
  # Add congruent sentence
  log.long$voice_congruent[log.long.start.row:log.long.end.row] <- log$incong_produced[row] 
  
  if(!is.na(log$pos[row])){
    log.long$pos[log.long.start.row:log.long.end.row] <- log$pos[row]
    if(log$pos[row]=='a'){
      log.long$onset_aux[log.long.start.row] <- 1
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
    }
    if(log$pos[row]=='v'){
      log.long$onset_verb[log.long.start.row] <- 1
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
    }
    ## BEGIN ADDED 9/21/2022
    if(log$pos[row]=='d'){
      log.long$onset_determiner[log.long.start.row] <- 1
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
    }
    if(log$pos[row]=='n'){
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
      ### BEGIN ADDED 10/21/2022
      if(log$substitution_error[row] %in% c('error','correction')){
        log.long$substitution_error[log.long.start.row] <- log$substitution_error[row]
        log.long$correct_noun[log.long.start.row] <- log$correct_noun[row]
        log.long$wrong_noun[log.long.start.row] <- log$wrong_noun[row]
      }
      ### END ADDED 10/21/2022
    }
    if(log$pos[row]=='b'){
      log.long$onset_by[log.long.start.row] <- 1
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
    }
    if(log$pos[row]=='p'){
      log.long$onset_passive_be[log.long.start.row] <- 1
      log.long$noun1[log.long.start.row] <- log$first_noun[row]
      log.long$noun2[log.long.start.row] <- log$second_noun[row]
      log.long$verb[log.long.start.row] <- log$verb[row]
      log.long$verb_voice[log.long.start.row] <- log$verb_voice[row]
    }
    ## END ADDED 9/21/2022
  }
  if(!is.na(log$dp_or_np[row])){log.long$dp_or_np[log.long.start.row:log.long.end.row] <- log$dp_or_np[row]}
  if(!is.na(log$role[row])){log.long$syntactic_role[log.long.start.row:log.long.end.row] <- log$role[row]}
  
  if((row%%100)==0){message("Time to complete row ",row," of ",nrow(log)," (",round(100*row/nrow(log),2),"%): ",round((Sys.time()-start.time),2))}
  if(!is.na(log.long$word[440018])){if(log.long$word[440018]=="getting"){message("'getting' at row: ",row)}}
}
message("Total big loop duration: ",round(Sys.time()-start.time,2)) # 7 mins

# Replace periods in column names with underscores
colnames(log.long) <- gsub(".","_",colnames(log.long),fixed=TRUE)

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
log.long$pos <- as.factor(log.long$pos)
log.long$dp_or_np <- as.factor(log.long$dp_or_np)
log.long$syntactic_role <- as.factor(log.long$syntactic_role)
log.long$word <- as.factor(log.long$word)

# Get rid of blank noun columns (mistakes in naming)
log.long$noun1[which(log.long$noun1=='')] <- NA
log.long$noun2[which(log.long$noun2=='')] <- NA
log.long$verb_voice[which(log.long$verb_voice=='')] <- NA
log.long$verb[which(log.long$verb=='')] <- NA

# Factorize
log.long$noun1 <- as.factor(log.long$noun1)
log.long$noun2 <- as.factor(log.long$noun2)
log.long$verb <- as.factor(log.long$verb)
log.long$verb_voice <- as.factor(log.long$verb_voice)
### BEGIN ADDED 10/21/2022
log.long$substitution_error <- as.factor(log.long$substitution_error)
log.long$correct_noun <- as.factor(log.long$correct_noun)
log.long$wrong_noun <- as.factor(log.long$wrong_noun)
### END ADDED 10/21/2022

# Write
save.dir <- paste0(path,patient,'/log/final/')
if(! dir.exists(save.dir)){
  dir.create(save.dir)
}
write.csv(log.long,
          paste0(save.dir,patient,'_log_512_aligned_with_ECoG_with_words.csv'),
          row.names=FALSE,
          quote=FALSE)

