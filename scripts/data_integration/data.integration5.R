## data.integration5.R
## David M. Budden 
## 26/10/2014

## Perform analysis for "Modelling transcriptional regulation of
## methylated and bivalent promoters" section

## packages
require(preprocessCore)
require(ggplot2)
require(gridExtra)
require(grid)

## parameters
data.path <- '../../results/'
results.path <- '../../results/'
cell.lines <- c('hesc', 'gm12878', 'k562')

num.perms <- 10
granularity <- 100

## functions
source('../common/adj.rsq.R')
source('../common/arsinh.R')

## define standard error function
se <- function(x) {
  sqrt(var(x)/length(x))
}

## define evaluate by split function
evaluate.by.split <- function(reg.element, all.elements, split.data) {
  ## benchmark
  benchmark.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(all.elements, collapse='+'))) 
  benchmark.model <- lm(benchmark.equation, data=split.data)
  ## evaluate against rep1
  benchmark.cor <- cor(benchmark.model$fitted.values, split.data$rpkm.rep1) # Pearson's correlation
  benchmark.rsq <- adj.rsq(y = split.data$rpkm.rep1, 
                     y.hat = benchmark.model$fitted.values, num.vars = length(all.elements))   
  ##print(sprintf('Benchmark: adj. R^2 = %.4f', benchmark.rsq))   
  
  ## order data by regulatory element of interest
  split.data <- split.data[with(split.data, order(split.data[,reg.element])),]
  
  ## containers  
  pos.rsq.vec <- vector(mode='numeric', length=granularity-2)
  neg.rsq.vec <- vector(mode='numeric', length=granularity-2)
  random.neg.vec <- vector(mode='numeric', length=granularity-2)
  random.pos.vec <- vector(mode='numeric', length=granularity-2)
  random.neg.se <- vector(mode='numeric', length=granularity-2)
  random.pos.se <- vector(mode='numeric', length=granularity-2)
  agr.adj.rsq <- vector(mode='numeric', length=granularity-2)
  
  ## store best split (+ve)
  best.pos.threshold <- 0; best.neg.threshold <- 0; best.agr.threshold <- 0
  best.pos.thres.rsq <- 0; best.neg.thres.rsq <- 0; best.agr.thres.rsq <- 0
  
  i <- 0
  for (threshold in floor(seq(from=1, to=nrow(split.data), 
                              length.out=granularity))[2:(granularity-1)]) {    
    i <- i+1
    
    ## split data based on threshold
    neg.data <- split.data[1:threshold,]
    pos.data <- split.data[(threshold+1):nrow(split.data),]
    
    ## train pos model
    pos.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(all.elements, collapse='+'))) 
    pos.model <- lm(pos.equation, data=pos.data)
    ## evaluate against rep1
    pos.cor <- cor(pos.model$fitted.values, pos.data$rpkm.rep1) # Pearson's correlation
    pos.rsq <- adj.rsq(y = pos.data$rpkm.rep1, 
                       y.hat = pos.model$fitted.values, num.vars = length(all.elements)) 
    
    ## train negative model
    neg.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(all.elements, collapse='+'))) 
    neg.model <- lm(neg.equation, data=neg.data)
    ## evaluate against rep1
    neg.cor <- cor(neg.model$fitted.values, neg.data$rpkm.rep1) # Pearson's correlation
    neg.rsq <- adj.rsq(y = neg.data$rpkm.rep1, 
                       y.hat = neg.model$fitted.values, num.vars = length(all.elements))
    
    ## aggregate performance
    agr.rsq <- adj.rsq(y=c(neg.data$rpkm.rep1, pos.data$rpkm.rep1),
              y.hat=c(neg.model$fitted.values, pos.model$fitted.values),
              num.vars = length(all.elements))
    
    ## train random models
    random.rsq <- vector(mode='numeric', length=num.perms)
    for (rand in 1:num.perms) {
      random.data <- split.data[sample(1:nrow(split.data), threshold, replace=TRUE),]
      
      random.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(all.elements, collapse='+'))) 
      random.model <- lm(random.equation, data=random.data)
      ## evaluate against rep1
      random.cor <- cor(random.model$fitted.values, random.data$rpkm.rep1) # Pearson's correlation
      random.rsq[rand] <- adj.rsq(y = random.data$rpkm.rep1, 
                                  y.hat = random.model$fitted.values, num.vars = length(all.elements))
    }    
    
    pos.rsq.vec[i] <- pos.rsq; neg.rsq.vec[i] <- neg.rsq; agr.adj.rsq[i] <- agr.rsq
    random.neg.vec[i] <- mean(random.rsq)
    random.pos.vec[granularity-i-1] <- mean(random.rsq)
    random.neg.se[i] <- se(random.rsq)
    random.pos.se[granularity-i-1] <- se(random.rsq)    
    
    ## DEBUG PRINT
    ##print(sprintf('pos(%d): %.4f; neg(%d): %.4f; (%.4f)', 
    ##              nrow(pos.data), pos.rsq, nrow(neg.data), neg.rsq, agr.adj.rsq[i]))
    
    ## check if best
    if (pos.rsq >= best.pos.thres.rsq) {
      best.pos.thres.rsq <- pos.rsq; best.pos.threshold <- threshold
    }
    if (neg.rsq >= best.neg.thres.rsq) {
      best.neg.thres.rsq <- neg.rsq; best.neg.threshold <- threshold
    }
    if (agr.rsq >= best.agr.thres.rsq) {
      best.agr.thres.rsq <- agr.rsq; best.agr.threshold <- threshold
    }
  }
  
  ## store results in single data frame
  results <- list()
  results[['summary']] <- data.frame(thresholds=floor(seq(from=1, to=nrow(split.data), 
                                                          length.out=granularity))[2:(granularity-1)],
                                      pos.rsq=pos.rsq.vec, neg.rsq=neg.rsq.vec, 
                                      rand.pos.rsq=random.pos.vec, rand.neg.rsq=random.neg.vec, 
                                      rand.pos.se=random.pos.se, rand.neg.se=random.neg.se,
                                      agr.rsq=agr.adj.rsq)
  results[['best.pos']] <- best.pos.threshold; results[['best.neg']] <- best.neg.threshold
  results[['best.agr']] <- best.agr.threshold; results[['benchmark']] <- benchmark.rsq
  
  return(results)
}

## MAIN CODE =====================================================================================>

for (cell.line in cell.lines) {
  ## prepare data
  histone.data <- read.csv(sprintf('%sintegration_1_%s.csv', data.path, cell.line))
  methyl.data <- read.csv(sprintf('%sintegration_3_%s.csv', data.path, cell.line))##[,
                    ##c('Ensembl.Gene.ID', 'rpkm.rep1', 'rpkm.rep2', 'MMFS.rep1', 'MMFS.rep2')]  
  prepared.data <- merge(methyl.data, histone.data[,c('Ensembl.Gene.ID', 'h2az.score', 'h3k4.score',
                    'h3k27.score', 'h3k9.score')], by='Ensembl.Gene.ID')
  all.elements <- c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score', 'MMFS.rep2')  
  
  
  ## FIRST SPLIT =================================================================================>    
  split1.data <- prepared.data
  split1.by <- 'MMFS.rep2'  
  split1.elements <- setdiff(all.elements, split1.by)
  #split1.elements <- all.elements
  
  ## run experiment
  results1 <- evaluate.by.split(split1.by, split1.elements, split1.data)
  
  ## sort (this happened inside the other function only)
  split1.data <- split1.data[with(split1.data, order(split1.data[,split1.by])),]
  split1.neg.data <- split1.data[1:results1$best.pos,]
  split1.pos.data <- split1.data[(results1$best.pos+1):nrow(split1.data),] 
  
  ## GENERATE PLOT FOR FIRST SPLIT ===============================================================>
  plot.data <- data.frame(pos.mean=results1$summary$pos.rsq-results1$benchmark, 
                          pos.se=results1$summary$rand.pos.se,
                          neg.mean=results1$summary$neg.rsq-results1$benchmark, 
                          neg.se=results1$summary$rand.neg.se,
                          agr.rsq=results1$summary$agr.rsq-results1$benchmark)
  
  split1.plot <- ggplot(plot.data, aes(x=1:nrow(plot.data), y=pos.mean)) + 
    geom_line(colour=rgb(0,0,1), size=1.0) +
    geom_errorbar(aes(ymin=pos.mean-pos.se, ymax=pos.mean+pos.se), colour=rgb(0,0,1)) +
    #geom_line(aes(y=neg.mean), colour=rgb(1,0,0), size=1.0) +
    #geom_errorbar(aes(ymin=neg.mean-neg.se, ymax=neg.mean+neg.se), colour=rgb(1,0,0)) +
    geom_hline(aes(yintercept=0), size=1.0) +
    #geom_vline(aes(xintercept=results1$best.neg/nrow(split1.data)*(granularity-1)), 
    #           color=rgb(1,0,0)) +
    geom_vline(aes(xintercept=results1$best.pos/nrow(split1.data)*(granularity-1)), 
               color=rgb(0,0,1), size=2) +
    #geom_vline(aes(xintercept=results1$best.agr/nrow(split1.data)*(granularity-1)), 
    #           color=rgb(0,0,0)) +
    #geom_line(aes(y=agr.rsq), colour=rgb(0,0,0), size=1.0) +
    ylab('Delta adj. R^2') +
    xlab('% total genes') +
    theme_bw(base_size=14)
  
  plot(split1.plot)
  ggsave(split1.plot, file=sprintf('%sfigure5_split1_%s.pdf', results.path, cell.line))
  
  ## SECONDS SPLIT ===============================================================================>      
  ## GET SECOND SPLIT AS NEGATIVE FROM FIRST 
  
  split2.data <- split1.neg.data
  split2.by <- 'h2az.score'  
  split2.elements <- setdiff(split1.elements, split2.by)
  #split2.elements <- split1.elements
  
  ## run experiment
  results2 <- evaluate.by.split(split2.by, split2.elements, split2.data)
  
  ## pick split 2 to maximise cumulative performance
  split2.data <- split2.data[with(split2.data, order(split2.data[,split2.by])),]
  split2.neg.data <- split2.data[1:results2$best.agr,]
  split2.pos.data <- split2.data[(results2$best.agr+1):nrow(split2.data),]
  
  ## GENERATE PLOT FOR SECOND SPLIT ==============================================================>
  plot2.data <- data.frame(pos.mean=results2$summary$pos.rsq-results2$benchmark, 
                          pos.se=results2$summary$rand.pos.se,
                          neg.mean=results2$summary$neg.rsq-results2$benchmark, 
                          neg.se=results2$summary$rand.neg.se,
                          agr.rsq=results2$summary$agr.rsq-results2$benchmark)
  
  split2.plot <- ggplot(plot2.data, aes(x=1:nrow(plot2.data), y=pos.mean)) + 
    geom_line(colour=rgb(0,0,1), size=1.0) +
    geom_errorbar(aes(ymin=pos.mean-pos.se, ymax=pos.mean+pos.se), colour=rgb(0,0,1)) +
    geom_line(aes(y=neg.mean), colour=rgb(1,0,0), size=1.0) +
    geom_errorbar(aes(ymin=neg.mean-neg.se, ymax=neg.mean+neg.se), colour=rgb(1,0,0)) +
    geom_hline(aes(yintercept=0)) +
    #geom_vline(aes(xintercept=results2$best.neg/nrow(split2.data)*(granularity-1)), 
    #           color=rgb(1,0,0)) +
    #geom_vline(aes(xintercept=results2$best.pos/nrow(split2.data)*(granularity-1)), 
    #           color=rgb(0,0,1)) +
    geom_vline(aes(xintercept=results2$best.agr/nrow(split2.data)*(granularity-1)), 
               color=rgb(0,0,0), size=2) +
    geom_line(aes(y=agr.rsq), colour=rgb(0,0,0), size=1.0) +
    ylab('Delta adj. R^2') +
    xlab('% total genes') +
    theme_bw(base_size=14)
  
  plot(split2.plot)
  ggsave(split2.plot, file=sprintf('%sfigure5_split2_%s.pdf', results.path, cell.line))
  
  ## FINAL SUMMARY (manuscript tree structure-specific)
  
  ## methylated genes
  meth.data <- split1.pos.data
  meth.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(split1.elements, collapse='+'))) 
  meth.model <- lm(meth.equation, data=meth.data)
  meth.rsq <- adj.rsq(y = meth.data$rpkm.rep1, 
              y.hat = meth.model$fitted.values, num.vars = length(split1.elements)) 
  
  ## bivalent genes
  bivalent.data <- split2.pos.data
  bivalent.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(split2.elements, collapse='+'))) 
  bivalent.model <- lm(bivalent.equation, data=bivalent.data)
  bivalent.rsq <- adj.rsq(y = bivalent.data$rpkm.rep1, 
              y.hat = bivalent.model$fitted.values, num.vars = length(split2.elements))
  
  ## non-bivalent genes
  active.data <- split2.neg.data
  active.equation <- as.formula(paste('rpkm.rep2 ~ ', paste(split2.elements, collapse='+'))) 
  active.model <- lm(active.equation, data=active.data)
  active.rsq <- adj.rsq(y = active.data$rpkm.rep1, 
                     y.hat = active.model$fitted.values, num.vars = length(split2.elements))
  
  ## write gene lists to file
  write.csv(meth.data, file=sprintf('%smmfs_pos_%s.csv', results.path, cell.line), row.names=F)
  write.csv(bivalent.data, file=sprintf('%sh2az_pos_%s.csv', results.path, cell.line), row.names=F)
  write.csv(active.data, file=sprintf('%sh2az_neg_%s.csv', results.path, cell.line), row.names=F)
  
  ## how did it go overall?
  
  print(sprintf('MMFS+ (%.0f%%): %.2f', 
                nrow(meth.data)/nrow(prepared.data)*100, meth.rsq-results1$benchmark))
  print(sprintf('H2A.Z+ (%.0f%%): %.2f', 
                nrow(bivalent.data)/nrow(prepared.data)*100, bivalent.rsq-results2$benchmark))
  print(sprintf('H2A.Z- (%.0f%%): %.2f', 
                nrow(active.data)/nrow(prepared.data)*100, active.rsq-results2$benchmark))
  
  write.csv(rbind(sprintf('MMFS+ (%.0f%%): %.2f', 
        nrow(meth.data)/nrow(prepared.data)*100, meth.rsq-results1$benchmark), 
        rbind(sprintf('H2A.Z+ (%.0f%%): %.2f', nrow(bivalent.data)/nrow(prepared.data)*100, 
        bivalent.rsq-results2$benchmark), sprintf('H2A.Z- (%.0f%%): %.2f', 
        nrow(active.data)/nrow(prepared.data)*100, active.rsq-results2$benchmark))),
        file=sprintf('%sTREEOME_results_%s.csv', results.path, cell.line), row.names=FALSE)
}
