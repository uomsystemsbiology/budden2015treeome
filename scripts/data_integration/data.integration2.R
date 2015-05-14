## data.integration2.R
## David M. Budden 
## 23/10/2014

## Perform analysis for 'Standard predictive modelling is unable to
## derive the regulatory signature of the H2A.Z histone variant' section

## packages
require(ggplot2)
require(gridExtra)

## parameters
data.path <- '../../results/'
results.path <- '../../results/'
cell.lines <- c('hesc', 'gm12878', 'k562')

## containers
results <- list()
adj.rsq <- list()
plots <- list()

## PERFORM PRINCIPAL COMPONENT ANALYSIS ==========================================================>

for (cell.line in cell.lines) {
  processed.data <- read.csv(sprintf('%sintegration_1_%s.csv', data.path, cell.line))
  hm.data <- processed.data[, c('h2az.score', 'h3k4.score', 'h3k27.score', 'h3k9.score')]
  colnames(hm.data) <- c('H2A.Z', 'H3K4me3', 'H3K27me3', 'H3K9me3')
  
  ## perform singular value decomposition
  svd.data <- svd(hm.data)
  u.matrix <- svd.data$u
  sigma.matrix <- diag(svd.data$d)
  v.matrix <- svd.data$v

  ## construct eigengene data
  pca.data <- data.frame(u.matrix %*% sigma.matrix)
  colnames(pca.data) <- paste('PC', c(1:ncol(pca.data)), sep='')
  
  ## evaluate predictive power of each PC
  pca.rsq <- vector()
  pca.t <- vector()
  for (i in 1:ncol(pca.data)) {
    model.rep1 <- lm(processed.data$rpkm.rep1 ~ pca.data[,i])
    pca.rsq[i] <- summary(model.rep1)$adj.r.squared
    pca.t[i] <- summary(model.rep1)$coefficients[2,3]
  }
  
  ## weighted loadings
  weighted.loadings <- v.matrix*sign(pca.t)
  results.summary <- data.frame(cbind(colnames(hm.data), weighted.loadings))
  colnames(results.summary) <- c('feature', paste('PC', c(1:ncol(pca.data)), sep=''))
  
  results[[cell.line]] <- results.summary
  adj.rsq[[cell.line]] <- pca.rsq
  
  
  ## GENERATE RESULTS FIGURES ====================================================================>
  
  ## relative loading of each HM to the PC most predictive of gene expression  
  features <- factor(colnames(hm.data), colnames(hm.data))  
  
  bar.plot <- ggplot(results.summary, aes(x=features, 
            y=as.numeric(as.character(results.summary[,sprintf('PC%s', which.max(pca.rsq))])))) +
    geom_bar(alpha=0.2, color=rgb(0,0,0), width=0.8, stat='identity') +
    ylab('weighted loading') +
    scale_y_continuous(limits=c(-1, 1)) +
    theme_bw(base_size=14) +
    theme(axis.text.x=element_text(angle=90, size=16, vjust=0.5, hjust=1), axis.title.x=element_blank()) +
    theme(axis.text.y=element_text(size=14), axis.title.y=element_text(size=20))  
  
  plot(bar.plot)
  ggsave(bar.plot, file=sprintf('%sfigure4_pca_%s.pdf', results.path, cell.line), width=6, height=6)
}

