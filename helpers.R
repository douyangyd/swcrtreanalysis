## Helper function: calculate CI with t-distribution
cal_confint_t <- function(pe, df, se){
  cihigh <- pe + qt(0.975, df)*se
  cilow <- pe - qt(0.975, df)*se
  return(cbind(cilow, cihigh))
}
