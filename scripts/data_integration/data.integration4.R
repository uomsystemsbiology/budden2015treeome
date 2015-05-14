## data.integration4.R
## David M. Budden 
## 25/10/2014

## Perform analysis for "naive predictive model integration is unsuitable for DNA 
## methylation data" section

## packages
require(preprocessCore)
require(ggplot2)
require(gridExtra)

## parameters
data.path <- '../../results/'
results.path <- '../../results/'
cell.lines <- c('hesc', 'gm12878', 'k562')

## functions
source('../common/adj.rsq.R')
source('../common/arsinh.R')

## containers
delta.rsq <- list()

for (cell.line in cell.lines) {
  ## prepare data
  histone.data <- read.csv(sprintf('%sintegration_1_%s.csv', data.path, cell.line))
  methyl.data <- read.csv(sprintf('%sintegration_3_%s.csv', data.path, cell.line))  
  prepared.data <- merge(histone.data[,c('Ensembl.Gene.ID', 'h2az.score', 'h3k4.score',
                  'h3k27.score', 'h3k9.score')], methyl.data, by='Ensembl.Gene.ID')
  
  ## build models
  
  ## evaluate performance before adding methylation scores  
  baseline <- as.formula(paste('rpkm.rep2 ~ ', paste(c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score'), collapse='+')))  
  model.baseline <- lm(baseline, data=prepared.data)
  ## evaluate against rep2
  cor.baseline <- cor(model.baseline$fitted.values, prepared.data$rpkm.rep1) # Pearson's correlation
  adj.rsq.baseline <- adj.rsq(y = prepared.data$rpkm.rep1, y.hat = model.baseline$fitted.values, num.vars = 4)
  
  ## container
  methyl.cor <- list()
  methyl.adj.rsq <- list()
  
  ## evaluate performance after adding methylation scores
  for (methyl.score in c('SMFS', 'MMFS', 'MMFR', 'SMRR')) {
    equation <- as.formula(paste('rpkm.rep2 ~ ', paste(c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score',
                                                 sprintf('%s.rep1', methyl.score)), collapse='+')))
    model <-  lm(equation, data=prepared.data)
    methyl.cor[[methyl.score]] <- cor(model$fitted.values, prepared.data$rpkm.rep1) # Pearson's correlation
    methyl.adj.rsq[[methyl.score]] <- adj.rsq(y = prepared.data$rpkm.rep1, y.hat = model$fitted.values, num.vars = 5)      
  }
  
  delta.rsq[[cell.line]] <- as.numeric(methyl.adj.rsq) - adj.rsq.baseline  
}

## here is how adding methylation data affects predictive power
print('Delta adj. R^2')

tmp <- rbind(delta.rsq$hesc, delta.rsq$gm12878)
results <- rbind(tmp, delta.rsq$k562)
print(results)

write.csv(results, file=sprintf('%sold_models_plus_meth.csv', results.path), row.names=c('h1-hESC', 'GM12878', 'K562'))

