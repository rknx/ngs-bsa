cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Starting R engine.'); cat('\n')

########## Environment preparation ##########

# Loading libraries
library("ggplot2")
library("ggrepel")
library("reshape2")

########## Input data ##########

# Read arguments from bash
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# Input candidate SNPs and all SNPs list
snp = read.delim(paste0(args[2],".ratio.tmp"), header=T)
if (args[6]=="y") cand = read.delim(paste0(args[2],".orf_SNPs.txt"), header =T)[,c("chr","pos","gene")]

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Finished data import.'); cat('\n')

########## Plotting functions ##########

manhattanplot = function(df) {

pl= ggplot(df, aes(pos, fit)) + #Position versus fitted value plot
	geom_point(size=0.4) + #Smoothing
	facet_grid (.~ chr, scales = "free_x", space = "free_x") + #Seperate graph by chromosome
	geom_hline(yintercept=as.numeric(args[3]),size=0.3,linetype="dotted") + #cut-off line
	geom_hline(yintercept=c(0.67,-0.67),size=0.3,linetype="dashed") + #significance threshold
	theme(legend.position="none") + #Remove legend
	scale_x_continuous(breaks=scales::pretty_breaks(n=36/length(unique(df$chr))), labels=sciNotn) + #X-axis formatting
    scale_y_continuous(breaks=seq(-1,1,0.1)) + #Y-axis formatting
	labs(x="Position", y="Ratio") #Axis names

if (args[6]=="y") pl = pl + geom_point(aes(pos, fit), df[df$gene,], shape=19) + geom_text_repel(aes(pos, fit, label=gene), df[df$gene,])
chrom = ifelse(length(unique(df$chr))==1, paste0("_",df$chr[1]), "")

ggsave(filename=paste0(args[2],chrom,"_loess_",args[3],"_",args[5],".pdf"), plot=pl, width=32, height=10, units="in")
ggsave(filename=paste0(args[2],chrom,"_loess_",args[3],"_",args[5],".png"), plot=pl, width=32, height=10, units="in")

}

allelegrid = function(df) {

pl= ggplot(df, aes(pos, value, col=bulk)) +
	geom_point(size=0.3, alpha=0.5) +
	facet_grid (.~ chr, scales = "free_x", space = "free_x") +
	scale_x_continuous(breaks=scales::pretty_breaks(n=36/length(unique(df$chr))), labels=sciNotn) +
	scale_y_continuous(breaks=seq(-1,1,0.1)) +
	labs(x="Position", y="Allele frequency")
	
if (args[6]=="y") pl=pl+geom_point(aes(x=pos, y=value), df[df$gene,], shape=19, stat = "unique")+geom_text_repel(aes(x=pos, y=value, label=gene), df[df$gene,], stat = "unique")
chrom = ifelse(length(unique(df$chr))==1, paste0("_",df$chr[1]), "")

ggsave(filename=paste0(args[2],chrom,"_allele_",args[3],"_",args[5],".pdf"), plot=pl, width=32, height=10, units="in")
ggsave(filename=paste0(args[2],chrom,"_allele_",args[3],"_",args[5],".png"), plot=pl, width=32, height=10, units="in")

}

########## Prettyplot functions ##########

# Scientific notation for x-axis
sciNotn = function(x, digits = 1) {
  print(x)
  expn = floor(log10(min(x[x>0])))
  base = round(x / 10^expn, digits)
  lab = vector()
  for(i in 1:length(x)) lab[i] = ifelse(x[i]==0,0,as.expression(substitute(base %*% 10^expn, list(base=base[i],expn=expn))))
  return (lab)
}

########## Data Preparation ##########

# Removing missing values
snp=snp[complete.cases(snp),]

# Spliting table by chromosomes
snps=split(snp, snp$chr)

snp=do.call(rbind,lapply(seq_along(snps), function(x) {
    df=snps[[x]]
    fit=loess(df$delratio~df$pos,degree=2,span=as.numeric(args[4]))$fitted
    absfit=loess(abs(df$delratio)~df$pos,degree=2,span=as.numeric(args[4]))$fitted
    return(data.frame(df,fit=fit,absfit=absfit))
  }))

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Finished LOESS smoothing.'); cat('\n')

# Merge candidates
if (args[6]=="y") snp=merge(snp,cand,all=T)

# Rename chromosomes
levels(snp[,1])=paste('Chromosome',1:length(unique(snp$chr)))

# save tables
save(snp,file=paste0(args[2],".snpdata.rda"))

# Prepare data for allele plot
if  (args[6]=="y") {
	snpm=melt(snp, id=c('chr','pos', 'delratio', 'fit', 'absfit', 'gene'), var='bulk')
} else {
	snpm=melt(snp, id=c('chr','pos', 'delratio', 'fit', 'absfit'), var='bulk')
}

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Ready to make plots.'); cat('\n')

######### Actual plotting ##########

# Getting all chromosome loess fitted plot
manhattanplot(snp)

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Manhattan plot generated for entire genome'); cat('\n')

# Plotting allele frequency
allelegrid(snpm)

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Allele ratio graph generated for entire genome.'); cat('\n')

# Plotting by chromosome
for (i in unique(snp$chr)) {
  manhattanplot(snp[snp$chr==i,]) # loess fitted plot
  allelegrid(snpm[snpm$chr==i,]) # allele frequency
}

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('All graphs generated for each chromosomes.'); cat('\n')

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'\t'); cat('Shutting down engine.'); cat('\n')

quit()
