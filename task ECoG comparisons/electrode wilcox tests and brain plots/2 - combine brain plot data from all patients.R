### Combine individual patient data into one big dataset
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
library('viridis') # for viridis color schemes
library('pals') # for ocean.balance, ocean.delta, ocean.curl colors
library('scico') # for palettes cork, vik, lisbon, tofino, berlin, roma


# Clean up
rm(list=ls())
cat("\014")

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 2
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 6
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 18
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 28
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


# Save path
save.data.path <- paste0(path,'analysis/R/task ECoG comparisons/electrode wilcox tests and brain plots/output/data/')

### Loop thru warped/unwarped data
for(warp.loop in c('warped','unwarped')){
  # warp.loop = c('warped','unwarped')[1]
  
  ### Loop thru frequency bands
  for(band.loop in c('high_gamma','beta')){
    # band.loop = c('high_gamma','beta')[1]
    
    ### Elec values
    for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
      # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[1]
      
      for(bin.loop in c(250,200,150,100,50)){ 
        # bin.loop <- 250
        
        data.dir <- paste0(save.data.path,
                           warp.loop,' data/',
                           band.loop,
                           '/patients - individually/',
                           lock.loop,'/',
                           bin.loop,'ms bins/')
        data.types <- list.files(data.dir)
        
        for(data.type.loop in data.types){ 
          # data.type.loop <- data.types[1]
          current.data.dir <- paste0(data.dir,data.type.loop,'/')
          current.files <- list.files(current.data.dir)
          
          all.data <- list()
          for(data.loop in 1:length(current.files)){ 
            # data.loop <- 1
            current.data <- current.files[data.loop]
            all.data[[data.loop]] <- read.csv(paste0(current.data.dir,current.data))
          }; rm(data.loop)
          
          # Combine
          all.data <- bind_rows(all.data)
          
          
          ### Save
          current.save.dir <- paste0(save.data.path,
                                     warp.loop,' data/',
                                     band.loop,
                                     '/patients - combined/',
                                     lock.loop,'/',
                                     bin.loop,'ms bins/')
          dir.create(current.save.dir, showWarnings = FALSE, recursive = TRUE)
          # Write data
          write.csv(all.data,
                    paste0(current.save.dir,data.type.loop,'.csv'),
                    row.names = FALSE, quote = FALSE)
          rm(all.data)
          
        }; rm(data.type.loop)
        message(warp.loop,' data - ',band.loop, ' - bin loop ',bin.loop,' - ',lock.loop, ' - complete! ',Sys.time())
      }; rm(bin.loop)
    }; rm(lock.loop)
  } # band.loop
} # warp.loop

message('Script completed successfully. ',Sys.time())
beep()

