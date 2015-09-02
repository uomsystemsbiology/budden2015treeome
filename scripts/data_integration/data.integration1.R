## data.integration1.R
## David M. Budden 
## 22/10/2014

## Perform analysis for 'Standard predictive modelling is unable to
## derive the regulatory signature of the H2A.Z histone variant' section

## packages
require(preprocessCore)
require(ggplot2)
require(gridExtra)
require(grid)

## parameters
data.path <- '../../data/'
results.path <- '../../results/'
cell.lines <- c('hesc', 'gm12878', 'k562')

## functions
source('../common/adj.rsq.R')
source('../common/arsinh.R')

## containers for storing useful output
scatter.plots <- list()
density.plots <- list()
processed.data <- list()

## repeat analysis for all cell lines of interest
for (cell.line in cell.lines) {
  ## data loci
  hm.path <- sprintf('%sencode/%s/hms/', data.path, cell.line)
  expression.path <- sprintf('%sencode/%s/expression/', data.path, cell.line)
  
  ## PRE-PROCESS HISTONE MODIFICATION DATA =======================================================>
  
  ## read histone modification data
  h2az.data <- read.csv(sprintf('%s%s_h2az.csv', hm.path, cell.line))[,c('Gene', 'name', 'hm.scores')]
  h3k4.data <- read.csv(sprintf('%s%s_h3k4me3.csv', hm.path, cell.line))[,c('Gene', 'name', 'hm.scores')]
  h3k27.data <- read.csv(sprintf('%s%s_h3k27me3.csv', hm.path, cell.line))[,c('Gene', 'name', 'hm.scores')]
  h3k9.data <- read.csv(sprintf('%s%s_h3k9me3.csv', hm.path, cell.line))[,c('Gene', 'name', 'hm.scores')]
  
  colnames(h2az.data) <- c('Ensembl.Gene.ID', 'name', 'h2az.score')
  colnames(h3k4.data) <- c('Ensembl.Gene.ID', 'name', 'h3k4.score')
  colnames(h3k27.data) <- c('Ensembl.Gene.ID', 'name', 'h3k27.score')
  colnames(h3k9.data) <- c('Ensembl.Gene.ID', 'name', 'h3k9.score')
  
  ## merge
  list.hm.data <- list(h2az.data, h3k4.data, h3k27.data, h3k9.data)
  hm.data <- Reduce(function(...) merge(..., all=T, by=c('Ensembl.Gene.ID', 'name')), list.hm.data)
  hm.data <- hm.data[!duplicated(hm.data$Ensembl.Gene.ID),]
  
  ## arsinh-transform
  hm.data[,c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score')] <-
    arsinh(hm.data[,c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score')])
  
  
  ## PRE-PROCESS EXPRESSION DATA =================================================================>
  
  ## read expression data
  expression.data.rep1 <- read.csv(sprintf('%sprocessed_output_%s_rep1.csv', 
                                    expression.path, cell.line))[,c('Ensembl.Gene.ID', 'counts', 'rpkm')]
  expression.data.rep2 <- read.csv(sprintf('%sprocessed_output_%s_rep2.csv', 
                                    expression.path, cell.line))[,c('Ensembl.Gene.ID', 'counts', 'rpkm')]
  
  colnames(expression.data.rep1) <- c('Ensembl.Gene.ID', 'counts.rep1', 'rpkm.rep1')
  colnames(expression.data.rep2) <- c('Ensembl.Gene.ID', 'counts.rep2', 'rpkm.rep2')
    
  ## ensure expression data is sorted
  expression.data.rep1 <- expression.data.rep1[with(expression.data.rep1, order(Ensembl.Gene.ID)),]
  expression.data.rep2 <- expression.data.rep2[with(expression.data.rep2, order(Ensembl.Gene.ID)),]
  
  ## merge
  expression.data <- cbind(expression.data.rep1, expression.data.rep2[,c('counts.rep2', 'rpkm.rep2')])
  expression.data <- expression.data[!duplicated(expression.data$Ensembl.Gene.ID),]
  
  ## remove low-quality RNA-seq data
  expression.data <- expression.data[(expression.data$counts.rep1 > 1) & (expression.data$counts.rep2 > 1),]
  
  ## quantile normalise and arsinh-transform
  expression.data[,c('rpkm.rep1', 'rpkm.rep2')] <-
    normalize.quantiles(as.matrix(expression.data[,c('rpkm.rep1', 'rpkm.rep2')]))
  expression.data[,c('rpkm.rep1', 'rpkm.rep2')] <-
    arsinh(expression.data[,c('rpkm.rep1', 'rpkm.rep2')])
  
  
  ## INTEGRATE DATA ==============================================================================>  
  
  ## merge expression and histone modification data
  prepared.data <- merge(expression.data[,c('Ensembl.Gene.ID','rpkm.rep1', 'rpkm.rep2')], hm.data, by='Ensembl.Gene.ID')
  processed.data[[cell.line]] <- prepared.data
  
  ## EVALUATE FIRST MODEL
  
  ## train regression model on rep1
  equation.rep1 <- as.formula(paste('rpkm.rep1 ~ ', paste(c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score'), collapse='+')))  
  model.rep1 <- lm(equation.rep1, data=prepared.data)
  ## evaluate against rep2
  cor.rep1 <- cor(model.rep1$fitted.values, prepared.data$rpkm.rep2) # Pearson's correlation
  adj.rsq.1 <- adj.rsq(y = prepared.data$rpkm.rep2, y.hat = model.rep1$fitted.values, num.vars = 4)
  
  ## EVALUATE SECOND MODEL
  equation.rep2 <- as.formula(paste('rpkm.rep2 ~ ', paste(c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score'), collapse='+')))  
  model.rep2 <- lm(equation.rep2, data=prepared.data)
  ## evaluate against rep1
  cor.rep2 <- cor(model.rep2$fitted.values, prepared.data$rpkm.rep1) # Pearson's correlation
  adj.rsq.2 <- adj.rsq(y = prepared.data$rpkm.rep1, y.hat = model.rep2$fitted.values, num.vars = 4)
  
  ## GENERATE RESULTS FIGURES ====================================================================>
  
  ## prepare data
  plot.data <- data.frame(measured=prepared.data$rpkm.rep1, predicted=model.rep2$fitted.values)
  
  ## model performance label
  label <- paste(sprintf('adj.~R^2 == %.2f', adj.rsq.2))  
  
  ## note: Figures are for replicate 2-trained models (the variation is minimal)
  scatter.plots[[cell.line]] <- ggplot(plot.data, aes(x=measured, y=predicted)) +
    geom_point(alpha=0.05, color=rgb(0,0,0)) +
    geom_smooth(method=lm, color=rgb(1,0,0), size=1, se=FALSE) +
    ylab('predicted expression') +
    xlab('RNA-seq expression') +
    scale_x_continuous(limits=c(0,10)) +
    scale_y_continuous(limits=c(0,6)) +
    annotate('text', x=2.5, y=5.5, label=label, parse=TRUE, size=6) +
    theme_bw(base_size=14)
  
  ##plot(scatter.plots[[cell.line]])
  
  density.plots[[cell.line]] <- ggplot(plot.data, aes(x=measured, y=1000*(..count../sum(..count..)))) +
    geom_density(color=rgb(0,0,0), fill=rgb(0,0,0), alpha=0.1) +
    theme_bw(base_size=14) +
    ylab('') + xlab('') +
    theme(axis.text.x=element_text(colour=rgb(1,1,1,0)), axis.text.y=element_text(colour=rgb(1,1,1,0))) +
    scale_x_continuous(limits=c(0,10)) +
    scale_y_continuous(limits=c(0,6))
}

## merge subpots
fig1 <- arrangeGrob(density.plots[['hesc']], density.plots[['gm12878']], density.plots[['k562']],
                    scatter.plots[['hesc']], scatter.plots[['gm12878']], scatter.plots[['k562']],
                    ncol=3, nrow=2, heights=c(1,3))
grid.draw(fig1)

## removing class check from ggsave as a workaround.  This will break at some point in the future
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

ggsave(fig1, file=sprintf('%sfigure3_regression.pdf', results.path), width=12.5, height=5.76)
## export at 8.97" by 5.76" (default)

## export data
write.csv(processed.data[['hesc']], file=sprintf('%sintegration_1_hesc.csv', results.path), row.names=FALSE)
write.csv(processed.data[['gm12878']], file=sprintf('%sintegration_1_gm12878.csv', results.path), row.names=FALSE)
write.csv(processed.data[['k562']], file=sprintf('%sintegration_1_k562.csv', results.path), row.names=FALSE)

