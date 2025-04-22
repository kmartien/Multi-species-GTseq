library(tidyverse)
library(ggplot2)
library(plyr)
library(gdata)
library(data.table)
library(stringr)
rm(list = ls())
######################################################################

args <- commandArgs(trailingOnly = TRUE) # for introducing parameters in command line
# arguments to add to command line: 
ID<-args[1] # species ID 
sample<-args[2] # sample ID
genome<-args[3] # genome ("FULL" or "RM" for repeat-masked)
winsize<-as.numeric(args[4]) #window size used for heterozygosity (100 or 1000 kb (1MB))
# ID<-"testsp"
# sample<-"z1234"
# genome<-"full"
# winsize<-1000

date<-format(Sys.Date(), "%d.%m.%y")

winhetfiles=list.files(pattern="_summary.het")
# extract the scaffold numbers (based on order in scaffold list, selected 
# by array range in bash script, e.g., chr2-5, 7-22 -> --array=2-5,7-22
Scaffolds2<-substring(winhetfiles,regexpr("_summary",winhetfiles) - 2)
Scaffolds3<-gsub("_summary.het","",Scaffolds2)
winhetfiles <- unique(grep(paste(Scaffolds3,collapse="|"),winhetfiles,value=TRUE))

# remove NA's 
winhetfiles<-winhetfiles[!is.na(winhetfiles)]
files <- lapply(winhetfiles, read.table, header=FALSE, sep="")
for (i in seq(files)) assign(paste0("Chr",Scaffolds3[[i]]),files[[i]])

dflist<-as.vector(ls()[sapply(ls(), function(i) class(get(i))) == "data.frame"])

allhet<-rbindlist(mget(ls(pattern = "Chr\\d+")), idcol = TRUE)

######################################################################
x <- allhet
colnames(x) <- c("chrom", "homo", "het", "calls", "prop_het")
x$V6<-ID
x$pct_het<-x$prop_het*100

pct_het_filter<-1.0 # filter for excluding window with >x% heterozygosity
calls_filter<-100000 # filter for minimum number of called sites per window

# Filter windows to remove those with min. called sites, 
# and remove windows with excess heterozygosity
pdf(file=paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_Het_plots_ANGSD_",date,".pdf"))
x %>%
  filter(prop_het<pct_het_filter) %>%
  filter(calls>calls_filter) %>% # filters out windows with <x sites called
  ggplot(aes(x=pct_het,y=V6, fill=V6)) +
  geom_violin(size=0.2, trim = TRUE, colour="blue", fill="blue") +
  xlab("percent heterozygosity") +
  labs(title=paste0(winsize,"kb window heterozygosity")) +
  labs(fill="Sample") +
  ylab(NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  theme(axis.title.x = element_text(face = "bold", size=12))+
  theme(axis.text.x = element_text(face = "bold", size=11, colour="black"))+
  theme(axis.text.y = element_text(face = "bold", size =11, colour="black"))+
  theme(axis.title.y = element_text(face = "bold", size=12)) 


x %>%
  filter(prop_het<pct_het_filter) %>%
  filter(calls>calls_filter) %>%
  ggplot(aes(x=het,y=V6, fill=V6)) +
  geom_violin(size=0.2, trim = TRUE, colour="blue", fill="blue") +
  xlab(paste0("Count of heterozygotes/",winsize,"kb")) +
  labs(title=paste0(winsize,"kb window heterozygosity")) +
  labs(fill="Sample") +
  ylab(NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  theme(axis.title.x = element_text(face = "bold", size=12))+
  theme(axis.text.x = element_text(face = "bold", size=11, colour="black"))+
  theme(axis.text.y = element_text(face = "bold", size =11, colour="black"))+
  theme(axis.title.y = element_text(face = "bold", size=12))
dev.off()

############# summary stats
# calculate minimum and maximum heterozygosity from all intervals
# overall heterozygosity (hets/calls across all scaffolds)

# first remove windows with pct_het>0.1 and <100k calls per 1MB window

y <- x %>% filter(pct_het<pct_het_filter) %>% filter(calls>calls_filter)

Het_win<-round(mean(y$prop_het),6)
sd_het <- round(sd(y$prop_het),6)

minHet_win<-round(min(y$het, na.rm=T),0)
maxhet_count<-round(max(y$het),0)
maxcalls_ID<-which.max(y$het)
Calls_maxhet<-y[maxcalls_ID,4]
Max_Het_win<-max(y$prop_het)

# total number of heterozygous sites
tot_hets <- round(sum(y$het),0)

# count number of windows
count_win<-nrow(y)
# Count number of windows with heterozygosity = 0
count_nohet<-sum(y$het < 0.99) # unclear where to draw the line when the het count is a fraction...
# percent of windows with no hets
percent_nohet<-count_nohet/count_win
# count number of windows with heterozygosity < 2 het/Mb)
count_het_2_per_MB<-sum(y$prop_het*1000000 < 2.01)
# count number of windows with het < 1/100kb
count_het_1_per_100kb<-sum(y$prop_het*100000 < 1.01)


# print out het count summary
a=paste0("Min. heterozygote count/",winsize,"kb = ",minHet_win)
b=paste0("Max. heterozygote count/",winsize,"kb = ",maxhet_count)
c=paste0("Total number of heterozygote sites = ", tot_hets)
d=paste0("Number of windows with het(0) = ",count_nohet)
e=paste0("Proportion of windows with het(0) = ",percent_nohet)
d1=paste0("Number of windows = ",count_win)
f=paste0("Max heterozygosity  = ", Max_Het_win)
g=paste0("Average heterozygosity = ", Het_win)
h=paste0("Standard deviation of het/window = ", sd_het)
i=paste0("Number of windows with het(<2/MB) = ",count_het_2_per_MB)
j=paste0("Number of windows with het(<1/100kb) = ",count_het_1_per_100kb)
k=paste0("maximum percent heterozygosity filter = ",pct_het_filter)
l=paste0("called genotypes filter = ",calls_filter)

pdf(paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_het_summ_stats_ANGSD_",date,".pdf"))
plot(NA, xlim=c(0,13), ylim=c(0,13), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,3,a, pos=4)
text(1,7,b, pos=4)
text(1,8,f, pos=4)
text(1,1,c, pos=4)
text(1,5,d1, pos=4)
text(1,6,d, pos=4)
text(1,4,e, pos=4)
text(1,9,g, pos=4)
text(1,2,h, pos=4)
text(1,10,i, pos=4)
text(1,11,j, pos=4)
text(1,12,k, pos=4)
text(1,13,l, pos=4)
points(rep(1,13),1:13, pch=15)
dev.off()


###############

# Get boundaries of chromosomes for plotting
pos=as.numeric(rownames(unique(data.frame(y$chrom)[1])))
pos=append(pos,length(y$chrom))
numpos=NULL
for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

# Set plot colors (2 colors alternating between chromosomes)
mycols=NULL
for (i in (seq(1,length(numpos), by=2))){mycols[i]="#2171b5"}
for (i in (seq(2,length(numpos), by=2))){mycols[i]="#6baed6"}

# Barplot of heterozygosity in windows across the genome
pdf(paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_winhet_barplot_",date,".pdf"), width=8, height=4, pointsize=10)
par(mar=c(8,5,2,1))
b=barplot(1000*y$het/y$calls, ylim=c(0,10*pct_het_filter),
          border=mycols[as.factor(y$chrom)], col=mycols[as.factor(y$chrom)], 
          ylab="Het / kb", , cex.lab=1.25, main=paste(ID,"(Het=",Het_win*100,"%)"))
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(y$chrom)), las=3, 
     line=-.25, cex.axis=.8)
dev.off()

# Barplot of heterozygosity in windows across the genome, without main label, scaffold labels
pdf(paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_winhet_barplot_NoTitle_",date,".pdf"), width=8, height=4, pointsize=14)
par(mar=c(8,5,2,1))
b=barplot(1000*y$het/y$calls, ylim=c(0,10*pct_het_filter),
          border=mycols[as.factor(y$chrom)], col=mycols[as.factor(y$chrom)], 
          ylab="Het / kb", cex.lab=1.25)
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(y$chrom)), las=3, 
     line=-.25, cex.axis=0.8)
dev.off()

# barplot of individual chromosome windows.
dflist = mget(ls(pattern = "Chr\\d+"))

dflist <- Map(cbind, dflist, chrom = names(dflist))


pdf(paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_chrom_barplot_",date,".pdf"), width=8, height=7, pointsize=14)
par(mar=c(8,5,2,1))
for(i in 1:length(dflist))
{
  d1=as.data.frame(dflist[[i]])
  barplot(1000*d1$V2/d1$V3, ylim=range(0,1000*d1$V2/d1$V3),#ylim=c(0,3),
          ylab="Het / kb", cex.lab=1.25, main=paste0(unique(d1$chrom)))
}
dev.off()

# Histogram of per-window heterozygosity
pdf(paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_winhet_hist_",date,".pdf"), width=4, height=4, pointsize=14)
par(mar=c(5,5,2,1))
h=hist(1000*y$het/y$calls, breaks=seq(0,5, by=0.1), ylim=c(0,sum(1000*y$het/y$calls)),
       border="#2171b5", col="#2171b5", ylab="# of windows", xlab="Het / kb", cex.lab=1.5,
       main=paste(ID,"(Het=",Het_win*100,"%)"))
dev.off()
# If error message = "some 'x' not counted; maybe 'breaks' do not span range of 'x'",
# then change the breaks range to include the max value of 1000*y$het/y$calls.

#### save workspace
save.image(file=paste0(ID,"_",sample,"_",genome,"_",winsize,"kb_ROH_plots_",date,".Rdata"))
