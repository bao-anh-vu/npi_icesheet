## Sequential mean calculation
function(curr_mean, new_value, n) {
    # n is the number of values seen so far
  return(curr_mean + (new_value - curr_mean) / (n+1))
}