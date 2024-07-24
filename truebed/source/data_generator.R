## data generator for CNN

data_generator <- function(input_paths, output_paths) {

  index <- 1
  
  function() {
    while(TRUE) {
      if (index > length(file_paths)) {
        index <<- 1
      }

        input <- readRDS(input_paths[index])
        output <- readRDS(output_paths[index])

        index <<- index + 1
        yield(list(input, output))

    }
  }
}