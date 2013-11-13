#Initializing libraries  ######################################################

source("http://bioconductor.org/biocLite.R")
setwd("~/B/Stud/02/Data/CEL")
library('limma')
library('affy')
library('annotate')
library("hgu95av2.db")

N = 1000	#Number of bootstrap iterations
M = 50	#Number of bootstrap sample;	Should be '20'/'50'
ADJ_TYPE = 'fdr'#Adjustment method;	Should be 'bonferroni'/'fdr'
REP = TRUE	#Replacement On/Off;	Should be on if M==20; off if M==50
RESULT_FILE_1 <- paste('1.0Simple_expression,bootstrap=', sprintf("%d", N),',sample=', sprintf("%d", M),',adj_type=', ADJ_TYPE, '.txt', sep = '')
RESULT_FILE_2 <- paste('1.0Genewise_statistics,bootstrap=', sprintf("%d", N),',sample=', sprintf("%d", M),',adj_type=', ADJ_TYPE, '.txt', sep = '')

celfiles <- system("ls *.CEL",intern=TRUE)
label = grepl('normal', celfiles)
label[label==TRUE] <- 'normal'
label[label==FALSE] <- 'tumor'
tumor_celfiles <- celfiles[label == 'tumor']
normal_celfiles <- celfiles[label == 'normal']

#Preprocessing data  ##########################################################

celfiles <- system("ls *.CEL",intern=TRUE)
affybatch <- ReadAffy(filenames=celfiles)
myRMA <- rma(affybatch)

#Probes to genes matching  ####################################################

tmp_exp <- exprs(myRMA)
exp_mat  <- as.matrix(tmp_exp)
probes <- rownames(exp_mat)

symbol <- hgu95av2SYMBOL
mapped_probes <- mappedkeys(symbol)
symbol <- as.list(symbol[mapped_probes])

missing <- probes[!probes%in%names(symbol)]
names(missing) <- missing
genesymbol <- append(symbol,missing)
genesymbols <- as.vector(unlist(genesymbol[probes]))

exp_mat <- data.frame(probes, genesymbols, exp_mat)
symbols <- as.vector(unique(exp_mat[,2]))

KK <- lapply(symbols, function(x) {     
cat("processing ",x, "\n")
     o <- subset(exp_mat,exp_mat[,2]==x)[,3:ncol(exp_mat)] 
     a <- apply(o, 2, median)
})

set <- t(do.call('cbind', KK))
rownames(set) <- symbols


#Bootstrap in Differeintial expression  #######################################

sign_vector = vector()
sign_gene_vector <- vector()

for (i in 1:N){

	print(i)

	i_tumor_celfiles <- sample(tumor_celfiles, M, replace = REP)
	i_normal_celfiles <- sample(normal_celfiles, M, replace = REP)
	i_tumor_celfiles <- unique(i_tumor_celfiles)
	i_normal_celfiles <- unique(i_normal_celfiles)
	i_celfiles <- c(i_tumor_celfiles, i_normal_celfiles)

	i_label = grepl('normal', i_celfiles)
	i_label[i_label==TRUE] <- 'normal'
	i_label[i_label==FALSE] <- 'tumor'

	i_set <- set[, i_celfiles]
	i_target <- cbind(i_celfiles, i_label)
	i_design <- cbind(WT = 1, slope = (i_target[, 2] == 'normal'))

	i_fit <- lmFit(i_set, i_design)
	i_fit <- eBayes(i_fit)

	i_result <- topTable(i_fit, coef=2, number=length(symbols), adjust.method=ADJ_TYPE)

	i_sign_result_tmp <- as.matrix(i_result)
	i_sign_result <- i_sign_result_tmp[as.numeric(i_sign_result_tmp[, 5]) <= 0.05,]
	i_sign_genes <- rownames(i_sign_result)
	i_result <- data.frame(symbols,i_result)

	print(table(i_result[,6]<=0.05))
	sign_vector <- c(sign_vector, table(i_result[,6]<=0.05)[2])
	sign_gene_vector <- c(sign_gene_vector, i_sign_genes)

}


#Saving the results  ##########################################################

print(sign_vector)
setwd("~/B/Stud/03/Results")
write.table(as.matrix(sign_vector), RESULT_FILE_1, row.names = FALSE, col.names = FALSE)
write.table(as.matrix(sign_gene_vector), RESULT_FILE_2, row.names = FALSE, col.names = FALSE)



