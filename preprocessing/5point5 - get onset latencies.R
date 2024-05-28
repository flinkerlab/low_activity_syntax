### Chicken Syntax
### Add production data to 512 Hz log
### May 2020
### adam.milton.morgan@gmail.com

rm(list=ls())

#####################
### Manually edit ###
#####################
patient <- 'Patient010'

#########################################
### Auto - manually check stuff below ###
#########################################

### Set path
if(Sys.info()['sysname'] == 'Darwin'){
  path = '/Users/am4611/Dropbox/Research/ChickenSyntax/data/'} # Mac
if(Sys.info()['sysname'] == 'Linux'){
  path = '/home/adam/Dropbox/Research/ChickenSyntax/data/'} # Ubuntu

# Read data
data.dir <- paste0(path,patient,'/log/final/')
log.512 <- read.csv(paste0(data.dir,patient,'_log_512_aligned_with_ECoG_with_words.csv'))

# Create function to visualize data
task.plot <- function(data, start.sample=1, end.sample=nrow(data), fast=FALSE){
  # Each naming trial part as a tall grey bar
  with(data[start.sample:end.sample,], 
       plot(x=time,
            y=onset_task_name_character, 
            col='grey', 
            type='l',
            xlab="time (ms)",
            yaxt='n',
            ylab='on/off',
            main='Experiment Events'))
  # Each trial onset as a tall black bar
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_trial, col='black'))
  # Each scene description/listing trial part as a 2/3 grey bar
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=(onset_task_describe_scene + onset_task_list_characters) / 1.5, col='grey'))
  
  # Each reading trial part as a 1/2 grey bar
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=(onset_task_read_sentence + onset_task_read_list) / 2, col='grey'))
  
  # Each production onset of a character that character's color
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_dracula / 3, col='navyblue'))
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_nurse / 3, col='turquoise'))
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_dog / 3, col='orange2'))
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_chicken / 3, col='gold'))
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_ninja / 3, col='darkred'))
  with(data[start.sample:end.sample,], 
       lines(x=time,
             y=onset_frankenstein / 3, col='green'))
  
  if(! fast){
    # Each image on screen as a 1/10 grey bar
    for(i in start.sample:end.sample){
      with(data[i,], segments(x0=time,
                              y0=0,
                              y1=image_on/10, 
                              col=if("dog" %in% image_character){'orange2'
                              }else if('chicken' %in% image_character){'gold'
                              }else if('dracula' %in% image_character){'navyblue'
                              }else if('ninja' %in% image_character){'darkred'
                              }else if('frankenstein' %in% image_character){'green'
                              }else if('nurse' %in% image_character){'turquoise'
                              }else{'darkgrey'}))}
    
    for(i in start.sample:end.sample){
      with(data[i,], segments(x0=time,
                              y0=0,
                              y1=as.numeric(! is.na(word))/20, 
                              col=if("dog" %in% word){'orange'
                              }else if('chicken' %in% word){'yellow'
                              }else if('dracula' %in% word){'blue3'
                              }else if('ninja' %in% word){'red'
                              }else if('frankenstein' %in% word){'chartreuse3'
                              }else if('nurse' %in% word){'cyan'
                              }else{'darkgrey'}))}
    
    text(x=data$time[start.sample + ((end.sample - start.sample) / 30)],
         y=c(1,.9,2/3,.5,.3,.1,.05)-.02,
         labels=c('new trial','naming task','picture description','reading','production onset','stim onset','word duration'),
         col=c('black','darkgrey','darkgrey','darkgrey','darkred','darkred','red'),
         pos=4)
  }
  # For aesthetics
  abline(h=0, col='white', lwd=2)  
}
task.plot(log.512, 502000, 547800, fast=TRUE)
task.plot(log.512, 1098960, 1106182, fast=FALSE)




### Add production latency for naming trials
log.512$samples_from_stimulus_onset <- NA
# Segment the task into chunks delimited by "onset_task_name_character"
onset.task.name.character.rows <- which(log.512$onset_task_name_character == 1)
n.name.character.tasks <- sum(log.512$onset_task_name_character)
warning.count <- 0
for(i in 1:n.name.character.tasks){ 
  # i=156 # for troubleshooting
  # Create temporary dataset that is just one naming trial
  # (as delineated by )
  start.row <- onset.task.name.character.rows[i]
  task.name.character.word <- log.512$image_character[start.row]
  end.row <- if(i == n.name.character.tasks){
    nrow(log.512)}else{onset.task.name.character.rows[i + 1]}
  temp.dat <- log.512[start.row:end.row,]
  #head(temp.dat[which(! is.na(temp.dat$word))-2,])
  # Remove rows after the next onset of a picture naming, scene description, new trial, scene reading, list description, list reading task
  rm(end.row)
  end.row <- with(temp.dat,
                  which(onset_task_name_character==1 |
                          ##onset_trial==1 |
                          onset_task_describe_scene==1 |
                          onset_task_read_question==1 |
                          onset_task_read_sentence==1 |
                          onset_task_see_arrow==1 |
                          onset_task_read_list==1 |
                          onset_task_list_characters==1))[2] - 1
  if(i == n.name.character.tasks){end.row <- nrow(temp.dat)}
  if(is.na(end.row)){end.row <- 2}
  # If there's a noun production in that window, add time from stimulus onset
  # Shorten
    temp.dat <- temp.dat[1:end.row,]
    
    # Find onset of the name of the character in the image
    current.character <- gsub("ecologically_valid_","",as.character(temp.dat$image_character[1]))
    naming.latency.in.samples <- 
      which(temp.dat[,paste0('onset_',current.character)] == 1 & 
              temp.dat$syntactic_role =='none')[1]

  if(is.na(naming.latency.in.samples)){
    # Print warning if there's no noun production in this window
    warning.count <- warning.count + 1
    message('Warning: Missing "',task.name.character.word,'" in naming task number ',i,' (rows ',start.row,' to ',start.row + end.row,'; time ',log.512$time[start.row],' to ',log.512$time[start.row + end.row],'), possibly because a trial was omitted due to patient mistake.')  
  }
  # Add to big dataframe
  log.512$samples_from_stimulus_onset[start.row:(start.row + nrow(temp.dat) - 1)] <- 0:(nrow(temp.dat)-1)
}; message('Total warnings: ',warning.count)


# Clean up
#rm(start.row, end.row, temp.dat, current.character, naming.latency.in.samples, i, n.name.character.tasks, onset.task.name.character.rows)

######################
### Manually check ###
######################

# Check
prod.onset <- which((log.512$onset_chicken == 1 |
                           log.512$onset_dog == 1 |
                           log.512$onset_dracula == 1 |
                           log.512$onset_frankenstein == 1 |
                           log.512$onset_nurse == 1 |
                           log.512$onset_ninja == 1) &
                          log.512$syntactic_role=='none')
i <- 4
log.512[(prod.onset[i]-1):(prod.onset[i]+1),c('image_character','task','onset_task_name_character','word','phase','onset_chicken','onset_dog','onset_dracula','onset_frankenstein','onset_ninja','onset_nurse','syntactic_role','samples_from_stimulus_onset')]
prod.latency <- as.numeric(as.character(log.512$samples_from_stimulus_onset[prod.onset[i]]))
log.512[(prod.onset[i]-prod.latency-1):(prod.onset[i]-prod.latency+1),c('image_character','task','onset_task_name_character','word','phase','onset_chicken','onset_dog','onset_dracula','onset_frankenstein','onset_ninja','onset_nurse','syntactic_role','samples_from_stimulus_onset')]

############
### Auto ###
############

### Add production latency for everything with 2 nouns: picture descriptions, reading sentences, listing, reading lists
# Segment the task into chunks delimited by "onset_task_name_character"
onset.two.character.task.rows <- 
  with(log.512, which(onset_task_describe_scene == 1 |
                        onset_task_read_sentence == 1 |
                        onset_task_list_characters == 1 |
                        onset_task_read_list == 1))
n.two.character.tasks <- length(onset.two.character.task.rows)
error.check <- 0
warning.count <- 0
for(i in 1:n.two.character.tasks){ # i=(n.two.character.tasks -1) # for troubleshooting
  error.check <- error.check + 1
  # i=13 # check Patient003 trials 13 and 15 for troubleshooting bc no objects
  # Create temporary dataset that is just one naming trial
  # (as delineated by )
  start.row <- onset.two.character.task.rows[i]
  end.row <- if(i == n.two.character.tasks){
    nrow(log.512)}else{onset.two.character.task.rows[i + 1]}
  temp.dat <- log.512[start.row:end.row,]
  
  # Remove everything after the "read sentence" task
  rm(end.row)
  end.row <- with(temp.dat,
                  which(#onset_trial==1 |
                          onset_task_name_character==1 |
                          onset_task_describe_scene==1 |
                          onset_task_read_question==1 |
                          onset_task_read_sentence==1 |
                          onset_task_see_arrow==1 |
                          onset_task_list_characters==1 |
                          onset_task_read_list==1))[2] - 1
  if(is.na(end.row)){end.row <- 2}
  temp.dat <- temp.dat[1:end.row,]
  
  # Find onset of the name of the subject and object nouns
  if(error.check %% 2==1){
    current.onset.subject.row <- 
      with(temp.dat, 
           which((task=='produce_list' & syntactic_role=='1' |
                    task=='produce_sentence' & syntactic_role=='s' |
                    task %in% c('read_sentence','read_list')) &
                   (onset_chicken==1 |
                      onset_dog==1 |
                      onset_dracula==1 |
                      onset_frankenstein==1 |
                      onset_ninja==1 |
                      onset_nurse==1)))
    current.onset.object.row <- 
      with(temp.dat, 
           which((task=='produce_list' & syntactic_role=='2' |
                    task=='produce_sentence' & syntactic_role%in%c('do','bo') |
                    task %in% c('read_sentence','read_list')) &
                   (onset_chicken==1 |
                      onset_dog==1 |
                      onset_dracula==1 |
                      onset_frankenstein==1 |
                      onset_ninja==1 |
                      onset_nurse==1)))
    if(length(current.onset.subject.row)==0){
      # Print warning if any
      warning.count <- warning.count + 1
      message('Warning: First noun missing in the two-character task (scene description, sentence reading, list production, or list reading) number ',i,' (rows ',start.row,' to ',start.row + end.row,'; time ',log.512$time[start.row],' to ',log.512$time[start.row + end.row],'), possibly because a trial was omitted due to patient mistake.')  
    }
    if(length(current.onset.object.row)==0){
      # Print warning if any
      warning.count <- warning.count + 1
      message('Warning: 2nd noun missing in the two-character task (scene description, sentence reading, list production, or list reading) number ',i,' (rows ',start.row,' to ',start.row + end.row,'; time ',log.512$time[start.row],' to ',log.512$time[start.row + end.row],'), possibly because a trial was omitted due to patient mistake.')  
    }
  }

  # Add to big dataframe
  log.512$samples_from_stimulus_onset[start.row:(start.row + nrow(temp.dat) - 1)] <- 0:(nrow(temp.dat)-1)
}; message('Total warnings: ',warning.count)

######################
### Manually check ###
######################

# Check
prod.onset <- which((log.512$onset_chicken == 1 |
                           log.512$onset_dog == 1 |
                           log.512$onset_dracula == 1 |
                           log.512$onset_frankenstein == 1 |
                           log.512$onset_nurse == 1 |
                           log.512$onset_ninja == 1) &
                          log.512$syntactic_role=='bo')
as.numeric(as.character(log.512$samples_from_stimulus_onset[prod.onset]))
which(is.na(as.numeric(as.character(log.512$samples_from_stimulus_onset[prod.onset]))))
#prod.onset <- which((log.512$onset_chicken == 1 | log.512$onset_dog == 1 | log.512$onset_dracula == 1 | log.512$onset_frankenstein == 1 | log.512$onset_nurse == 1 | log.512$onset_ninja == 1) & log.512$phase == "character_naming_quick")
i <- 6
log.512[(prod.onset[i]-1):(prod.onset[i]+1),c('image_action','task','onset_task_describe_scene','onset_task_name_character','word','phase','onset_chicken','onset_dog','onset_dracula','onset_frankenstein','onset_ninja','onset_nurse','syntactic_role','samples_from_stimulus_onset')]
prod.latency <- as.numeric(as.character(log.512$samples_from_stimulus_onset[prod.onset[i]]))
log.512[(prod.onset[i]-prod.latency-1):(prod.onset[i]-prod.latency+1),c('image_action','task','onset_task_describe_scene','word','phase','onset_chicken','onset_dog','onset_dracula','onset_frankenstein','onset_ninja','onset_nurse','syntactic_role','samples_from_stimulus_onset')]

############
### Auto ###
############

### Write
# Write over the previous file
write.csv(log.512,
          paste0(data.dir,patient,'_log_512_aligned_with_ECoG_with_words_with_latencies.csv'),
          row.names=FALSE,
          quote=FALSE)


