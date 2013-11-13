#vect <- read.table()
#v <- as.matrix(vect)

dx <- 10

setwd("~/B/Stud/03/Results")
filename <- c('1.0Simple_expression,bootstrap=1000,sample=50,adj_type=fdr.txt','1.0Simple_expression,bootstrap=1000,sample=50,adj_type=bonferroni.txt','1.0Simple_expression,bootstrap=1000,sample=20,adj_type=fdr.txt','1.0Simple_expression,bootstrap=1000,sample=20,adj_type=bonferroni.txt')

for (f in filename){
	print(f)

	v <- read.table(f)
	v <- as.matrix(v)
	v <- v[!is.na(v)]
	u <- unique(v)

	brk <- seq(min(u)-dx, max(u)+dx,dx)
	postscript(paste(f,'.eps',sep=''),horizontal = FALSE, onefile = FALSE)
	hist(v, breaks = brk,main = f, col = 'red', border = 'red')
	dev.off()
}
warnings()
