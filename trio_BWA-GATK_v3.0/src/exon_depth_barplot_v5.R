args = commandArgs(trailingOnly=TRUE)

plot <- function(input_f,gene,sample){
	#x_name = paste("exon_",seq(nrow(input_f)),sep="")
	x_name = input_f$V1
	max = max(input_f$V2)+50
	p = barplot(input_f$V2, space=2, col="red", main=paste(gene,"exon region count of",sample), ylim=c(0,max))
	text(p, input_f$V2, labels = x_name, pos=3)
}

if (length(args)!=1){
        stop("missing the argument for gene name", call.=FALSE)
} else {
        gene = args[1]
}

read.table("id")->id
id=as.vector(unlist(id))

nsamples <- length(id)
pdf(paste0(gene,"_exon_depth.pdf"), width=10,height=10*nsamples)
#par(mfrow=c(ceiling(nsamples/2),2)) #ceiling takes a single numeric argument x and returns a numeric vector containing the smallest integers not less than the corresponding elements of x
#par(mfrow=c(nsamples,1), mar=c(1,1,1,1))
par(mfrow=c(nsamples,1))
for (i in 1:nsamples) {
	tab = read.table(paste0(id[i],"_",gene,"_exon.depth.tsv"))
	gene = args[1]
	sample = id[i]
	plot(tab,gene,sample)
}
dev.off()


