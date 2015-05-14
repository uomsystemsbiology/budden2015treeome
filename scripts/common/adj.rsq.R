## adj.rsq.R
## David M. Budden 
## 18/11/2014

adj.rsq <- function(y, y.hat, num.vars) {
  1-(1-cor(y, y.hat)^2)*((length(y)-1)/(length(y)-num.vars-1))
}
