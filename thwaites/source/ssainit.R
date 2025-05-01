ssainit <- function(p, x, ug, drivestress, gamma, initchoice) {
  if (missing(p) || missing(x) || missing(drivestress) || missing(gamma) || missing(initchoice)) {
    stop("All 5 input arguments are required.")
  }

  if (initchoice == 1) {
    # Ice shelves: linear profile from calving front condition
    u0 <- ug + gamma * x

  } else if (initchoice == 2) {
    # Ice streams: velocity depends on driving stress only
    u0 <- (-drivestress / p$C)^(1 / p$m)

  } else if (initchoice == 3) {
    # Predefined linear solution
    u0 <- (30.0 / p$secpera) * rep(1, length(x)) +
          (30.0 / p$secpera) * (x / p$L)

  } else {
    stop("initchoice must be 1, 2, or 3.")
  }

  return(u0)
}
