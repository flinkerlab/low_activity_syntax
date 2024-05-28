### Time convert function
### Winter 2021
### adam.milton.morgan@gmail.com


### Formats for input/output (e.g.):
# 'samples': -512 or 1024
# 'sample.labels': "sample_neg512" or "sample_1024"
# 'times': -1000 or 2000 

time.convert <- function(input, 
                         input.type = 'samples', 
                         output.type = 'times',
                         drop.labels = FALSE){
  if(input.type == 'sample.labels'){
    sample.labels = input
    samples = as.numeric(gsub('sample_','',gsub('neg','-',sample.labels)))
    times = round(samples / 512 * 1000)
  }
  if(input.type == 'samples'){
    samples = as.numeric(input)
    sample.labels = paste0('sample_',gsub('-','neg',as.character(samples)))
    times = round(samples / 512 * 1000)
  }
  if(input.type == 'times'){
    times = as.numeric(input)
    samples = round(times * 512 / 1000)
    sample.labels = paste0('sample_',gsub('-','neg',as.character(samples)))
  }
  output <- list('times'=times,
                 'samples'=samples, 
                 'sample.labels'=sample.labels)[[output.type]]
  if((! is.null(names(input))) & (! drop.labels)){
    names(output) <- names(input)
  }
  return(output)
}