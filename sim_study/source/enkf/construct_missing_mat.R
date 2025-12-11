## Construct missing observation matrix
construct_missing_matrix <- function(missing_pattern){
      
      missing_mat_ls <- lapply(1:ncol(missing_pattern), function(y) {
        mp_year <- missing_pattern[, y]
        nonmiss <- sum(mp_year)
        C_t <- sparseMatrix(i = 1:nonmiss, j = which(mp_year == 1), x = 1,
                            dims = c(nonmiss, nrow(missing_pattern)))
        return(C_t)
      })

      return(missing_mat_ls)
  }