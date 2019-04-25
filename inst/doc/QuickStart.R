## ----load-packages, message=FALSE, warning=F-----------------------------
require(BALLI)

## ------------------------------------------------------------------------
GenerateData <- function(nRow) {
	expr_mean <- runif(1,10,100)
	expr_size <- runif(1,1,10)
	expr <- rnbinom(20,mu=expr_mean,size=expr_size)
	return(expr)
}

data <- data.frame(t(sapply(1:10000,GenerateData)))
colnames(data) <- c(paste0("A",1:10),paste0("B",1:10))
rownames(data) <- paste0("gene",1:10000)
head(data)

## ------------------------------------------------------------------------
Group <- c(rep("A",10),rep("B",10))
Group

## ------------------------------------------------------------------------
design <- model.matrix(~Group, data = data)
head(design)

## ------------------------------------------------------------------------
dge <- DGEList(counts=data, group=Group)
dge <- calcNormFactors(dge)
dge

## ------------------------------------------------------------------------
tV <- tecVarEstim(dge,design)
tV

## ------------------------------------------------------------------------
fit <- balli(tV,intV=2)
fit

