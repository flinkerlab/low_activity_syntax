### Average across dataframes/matrices by element 
### Spring 2022
### adam.milton.morgan@gmail.com

# For 3D data (time x electrode x trial), a wrapper to perform some function across one of those dimensions (e.g., to get the mean trial activity for each electrode)

elementwise.matrix.apply <- function(list.of.matrices,
                                     .function = c("mean","sd")[1]){
  
  # Make sure matrices are data.frames
  if("matrix" %in% unlist(sapply(list.of.matrices, class))){
    list.of.matrices <- lapply(list.of.matrices, function(x){
      x <- data.frame(x)
      return(x)
    })
  } # If not a dataframe
  
  # Get dimensions
  n.rows <- nrow(list.of.matrices[[1]])
  n.cols <- ncol(list.of.matrices[[1]])
  
  # Get function
  .function <- get(.function)

  # Initialize output
  output <- data.frame(matrix(nrow = n.rows, ncol = n.cols))
  colnames(output) <- colnames(list.of.matrices[[1]])
  rownames(output) <- rownames(list.of.matrices[[1]])
  
  # Convert to 3D arrayy
  library('str2str')
  list.of.matrices <- ld2a(list.of.matrices) 
  
  # # Loop thru cells and apply function across data.frames
  # for(row.loop in 1:n.rows){ # row.loop = 1
  #   for(col.loop in 1:n.cols){ # col.loop = 2
  #     output[row.loop, col.loop] <-
  #       .function(list.of.matrices[row.loop, col.loop, ])
  #   } # col.loop
  # } # row.loop
  
  # Faster:
  output <- apply(list.of.matrices, c(1,2), .function)
  
  # Output
  return(output)
}  