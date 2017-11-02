errorBarSeg <- function(x,y,sd) {
  plot (x, y)
  segments(x,y-sd,x,y+sd)
  epsilon <- 0.02
  segments(x-epsilon,y-sd,x+epsilon,y-sd)
  segments(x-epsilon,y+sd,x+epsilon,y+sd)
}
