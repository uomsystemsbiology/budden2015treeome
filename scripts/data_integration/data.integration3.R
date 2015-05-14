## data.integration3.R
## David M. Budden 
## 25/10/2014

## Perform analysis for 'MMFS-quantified promoter methylation is 
## anti-correlated with gene expression' section

## packages
require(preprocessCore)
require(ggplot2)
require(gridExtra)

## parameters
data.path <- '../../data/'
results.path <- '../../results/'
cell.lines <- c('hesc', 'gm12878', 'k562')

## functions
source('../common/adj.rsq.R')
source('../common/arsinh.R')

## containers
processed.data <- list()

for (cell.line in cell.lines) {
  
  ## PRE-PROCESS EXPRESSION DATA =================================================================>
  
  ## read expression data
  expression.path <- sprintf('%sencode/%s/expression/', data.path, cell.line)
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
  
  
  ## PRE-PROCESS METHYLATION DATA ================================================================>    
  
  ## read methylation data
  methylation.path <- sprintf('%s/encode/%s/methylation/', data.path, cell.line)
  methylation.rep1 <- read.csv(sprintf('%s%s_methylation_scores_replicate_1_2000bp.csv',
                        methylation.path, cell.line, cell.line))[,c('ensemblID', 'OGS',
                        'sumMethylationFractionBySite', 'meanMethylationFractionBySite',
                        'meanMethylationFractionByRegion', 'sumOfScaledMethylatedReadsByRegion')]
  methylation.rep2 <- read.csv(sprintf('%s%s_methylation_scores_replicate_2_2000bp.csv',
                        methylation.path, cell.line, cell.line))[,c('ensemblID', 'OGS',
                        'sumMethylationFractionBySite', 'meanMethylationFractionBySite',
                        'meanMethylationFractionByRegion', 'sumOfScaledMethylatedReadsByRegion')]  
  colnames(methylation.rep1) <- c('Ensembl.Gene.ID', 'name', 'SMFS.rep1', 'MMFS.rep1', 'MMFR.rep1', 'SMRR.rep1')
  colnames(methylation.rep2) <- c('Ensembl.Gene.ID', 'name', 'SMFS.rep2', 'MMFS.rep2', 'MMFR.rep2', 'SMRR.rep2')
  
  ## ensure methylation data is sorted
  methylation.rep1 <- methylation.rep1[with(methylation.rep1, order(Ensembl.Gene.ID)),]
  methylation.rep2 <- methylation.rep2[with(methylation.rep2, order(Ensembl.Gene.ID)),]
  
  ## merge
  methylation.data <- cbind(methylation.rep1, methylation.rep2[,c('SMFS.rep2', 'MMFS.rep2', 'MMFR.rep2', 'SMRR.rep2')])
  methylation.data <- methylation.data[!duplicated(methylation.data$Ensembl.Gene.ID),]
  
  ## remove low-quality data
  methylation.data <- methylation.data[(methylation.data$SMFS.rep1 > 0) & (methylation.data$SMFS.rep2 > 0),]
  
  ## pairwise quantile normalise and arsinh-transform
  methylation.data[,c('SMFS.rep1', 'SMFS.rep2')] <-
    normalize.quantiles(as.matrix(methylation.data[,c('SMFS.rep1', 'SMFS.rep2')]))
  methylation.data[,c('MMFS.rep1', 'MMFS.rep2')] <-
    normalize.quantiles(as.matrix(methylation.data[,c('MMFS.rep1', 'MMFS.rep2')]))
  methylation.data[,c('MMFR.rep1', 'MMFR.rep2')] <-
    normalize.quantiles(as.matrix(methylation.data[,c('MMFR.rep1', 'MMFR.rep2')]))
  methylation.data[,c('SMRR.rep1', 'SMRR.rep2')] <-
    normalize.quantiles(as.matrix(methylation.data[,c('SMRR.rep1', 'SMRR.rep2')]))
  
  methylation.data[,c('SMFS.rep1', 'SMFS.rep2', 'MMFS.rep1', 'MMFS.rep2', 
                      'MMFR.rep1', 'MMFR.rep2', 'SMRR.rep1', 'SMRR.rep2')] <-
    arsinh(methylation.data[,c('SMFS.rep1', 'SMFS.rep2', 'MMFS.rep1', 'MMFS.rep2', 
                               'MMFR.rep1', 'MMFR.rep2', 'SMRR.rep1', 'SMRR.rep2')])
  
  ## DATA INTEGRATION ============================================================================> 
  
  ## merge data
  processed.data[[cell.line]] <- merge(expression.data, methylation.data, by='Ensembl.Gene.ID')
  
  ## write to file
  write.csv(processed.data[[cell.line]], file=sprintf('%sintegration_3_%s.csv', results.path, cell.line), row.names=F)
}

## analyse methylation/expression correlations

hesc.results <- as.data.frame(cor(processed.data[['hesc']][,c('rpkm.rep1', 'rpkm.rep2')], processed.data[['hesc']][,
                c('SMFS.rep1','SMFS.rep2','MMFS.rep1','MMFS.rep2','MMFR.rep1','MMFR.rep2','SMRR.rep1','SMRR.rep2')]))
gm12878.results <- as.data.frame(cor(processed.data[['gm12878']][,c('rpkm.rep1', 'rpkm.rep2')], 
                processed.data[['gm12878']][,c('SMFS.rep1','SMFS.rep2','MMFS.rep1','MMFS.rep2',
                'MMFR.rep1','MMFR.rep2','SMRR.rep1','SMRR.rep2')]))
k562.results <- as.data.frame(cor(processed.data[['k562']][,c('rpkm.rep1', 'rpkm.rep2')], 
                processed.data[['k562']][,c('SMFS.rep1','SMFS.rep2','MMFS.rep1','MMFS.rep2',
                'MMFR.rep1','MMFR.rep2','SMRR.rep1','SMRR.rep2')]))

row.names(hesc.results) <- c('H1-hESC rep1', 'H1-hESC rep2')
row.names(gm12878.results) <- c('GM12878 rep1', 'GM12878 rep2')
row.names(k562.results) <- c('K562 rep1', 'K562 rep2')

results <- rbind(hesc.results, gm12878.results)
results <- rbind(results, k562.results)

## visualise correlaton for best performing score

plot.cell.line <- 'k562'
plot.methyl <- 'MMFS.rep2'

label <- paste(sprintf('adj.~R^2 == %.2f', 
          adj.rsq(processed.data[[plot.cell.line]]$rpkm.rep1, processed.data[[plot.cell.line]][,plot.methyl], 1)))

scatter.plot <- ggplot(processed.data[[plot.cell.line]], aes(x=rpkm.rep1, y=processed.data[[plot.cell.line]][,plot.methyl])) +
  geom_point(alpha=0.05, color=rgb(0,0,0)) +
  geom_smooth(method=lm, color=rgb(1,0,0), size=1, se=FALSE) +
  ylab('Promoter methylation') +
  xlab('RNA-seq expression') +
  scale_x_continuous(limits=c(0,10)) +
  scale_y_continuous(limits=c(0,6)) +
  annotate('text', x=7.5, y=5.75, label=label, parse=TRUE, size=6) +
  theme_bw(base_size=14)

plot(scatter.plot)
print(results)

write.csv(results, file=sprintf('%smethyl_score_comparison.csv', results.path))
