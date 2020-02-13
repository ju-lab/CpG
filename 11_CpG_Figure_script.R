##Figure 1e
setwd("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision")

#no differentiate exon and no-exon
fig1e <- read.table("figure1e_reclassfication.txt",header=T)
fig1e <- t(fig1e)
fig1e<- as.data.frame(fig1e)
colnames(fig1e) <- c("Total","Thousand","PCAWG","Union")
fig1e
fig1e$ratio_thousand <- fig1e$Thousand/fig1e$Total
fig1e$ratio_pcawg <- fig1e$PCAWG/fig1e$Total
fig1e$ratio_union <- fig1e$Union/fig1e$Total

fig1e_reaarange <- fig1e[c(5,2,4,6,3,1),]
par(mfrow=c(1,1))
pdf("Figure1e_union_errorbar.pdf")
barplot(fig1e_reaarange$ratio_thousand,ylim=c(0,0.30),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'Thousand' )
barplot(fig1e_reaarange$ratio_pcawg,ylim=c(0,0.30),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'PCAWG' )
barplot(fig1e_reaarange$ratio_union,ylim=c(0,0.30),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'Union' )
dev.off()

#add confidence interval using prop.test (alpha = 0.05)
fig1e_reaarange$lower <- 'NA'
fig1e_reaarange$upper <- 'NA'
for(i in 1:6){
  fig1e_reaarange[i,]$lower<-prop.test(fig1e_reaarange[i,]$Union,fig1e_reaarange[i,]$Total,conf.level = 0.95,correct=F)$conf.int[1]
  fig1e_reaarange[i,]$upper<-prop.test(fig1e_reaarange[i,]$Union,fig1e_reaarange[i,]$Total,conf.level = 0.95,correct=F)$conf.int[2]
}
fig1e_reaarange

pdf("Figure1e_union_errorbar.pdf")
barplot(fig1e_reaarange$ratio_union,ylim=c(0,0.30),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'Union' )
arrows(a[,1],as.numeric(fig1e_reaarange$lower),a[,1],as.numeric(fig1e_reaarange$upper),lwd=1.5,angle=90,code=3,length=0.05,col=c('red','orange','blue','skyblue','green','yellow'))
dev.off()

#number of samples in each cohort
temp<- read.table("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/PCAWG-8/temp.txt",header=F)
ncol(temp)
temp[1:5,1:10]
2504+1823


####coding and non-coding####
#exonic
fig1e_exonic <- read.table("figure1e_exonic_reclassfication.txt",header=T)
fig1e_exonic<- t(fig1e_exonic)
fig1e_exonic
colnames(fig1e_exonic) <- c("Total","Thousand","PCAWG","Union")
fig1e_exonic <- as.data.frame(fig1e_exonic)
fig1e_exonic$ratio_thousand <- fig1e_exonic$Thousand/fig1e_exonic$Total
fig1e_exonic$ratio_pcawg <- fig1e_exonic$PCAWG/fig1e_exonic$Total
fig1e_exonic$ratio_union <- fig1e_exonic$Union/fig1e_exonic$Total
fig1e_exonic_rearrange <- fig1e_exonic[c(6,7,2,5,8,9,4,1,3,10,11),]
fig1e_exonic_rearrange
pdf("fig1e_exonic_difference.pdf")
par(mfrow=c(1,1))
barplot(fig1e_exonic_rearrange$ratio_union,ylim=c(0,0.30),col=c('red','red4','orange','blue','skyblue','skyblue4','green','yellow','yellow4',"darkseagreen1","darkseagreen4"),names.arg=c('TSS_coding_noexon','TSS_coding','TSS_noncoding_noexon','Intergenic','Intragenic_coding_noexon','Intragenic_coding','Intragenic_noncoding','Not in CGI_noexon','Not in CGI','mis_noexon','mis_exon'),space=c(1,0.2,1,1,1,0.2,1,1,0.2,1,0.2),xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','red4','orange','blue','skyblue','skyblue4','green','yellow','yellow4',"darkseagreen1","darkseagreen4"),cex.names = 0.8, main = 'Union')
dev.off()

#syn,nonsyn,stop,unk,noexon
fig1e_exonic <- read.table("figure1e_reclassfication_selectivepressure.txt",header=T)
fig1e_exonic<- t(fig1e_exonic)
fig1e_exonic
colnames(fig1e_exonic) <- c("Total","Thousand","PCAWG","Union")
fig1e_exonic <- as.data.frame(fig1e_exonic)
fig1e_exonic$ratio_thousand <- fig1e_exonic$Thousand/fig1e_exonic$Total
fig1e_exonic$ratio_pcawg <- fig1e_exonic$PCAWG/fig1e_exonic$Total
fig1e_exonic$ratio_union <- fig1e_exonic$Union/fig1e_exonic$Total
fig1e_exonic_rearrange <- fig1e_exonic[c(7,9,8,10,2,6,12,13,14,15,5,1,3,4,11),]
fig1e_exonic_rearrange
pdf("fig1e_selective_pressure.pdf")
par(mfrow=c(1,1))
barplot(fig1e_exonic_rearrange$ratio_union,ylim=c(0,0.30),col=c('red1','red2','red3','red4',
                                                                'orange','blue','skyblue1','skyblue2','skyblue3','skyblue4',
                                                                'green','yellow1','yellow2','yellow3','yellow4'),
        names.arg=c('TSS_coding_noexon','TSS_coding_syn','TSS_coding_non','TSS_coding_sto',
                    'TSS_noncoding_noexon','Intergenic',
                    'Intragenic_coding_noexon','Intragenic_coding_syn','Intragenic_coding_non','Intragenic_coding_sto',
                    'Intragenic_noncoding',
                    'Not in CGI_noexon','Not in CGI_syn','Not in CGI_non','Not in CGI_sto'),
        space=c(1,0.2,0.2,0.2,1,1,1,0.2,0.2,0.2,1,1,0.2,0.2,0.2),xlab = 'Types of CpG islands', ylab='Mutation rate',
        border=c('red1','red2','red3','red4','orange','blue','skyblue1','skyblue2','skyblue3','skyblue4','green','yellow1','yellow2','yellow3','yellow4'),
        cex.names = 0.8, main = 'Union')
dev.off()


### Thousand genome race diffence??
fig1e <- read.table("figure1e_race.txt",header=T)
fig1e <- t(fig1e)
fig1e<- as.data.frame(fig1e)
colnames(fig1e) <- c("Total","EAS","AMR","AFR","EUR","SAS")
fig1e
fig1e$ratio_EAS <- fig1e$EAS/fig1e$Total
fig1e$ratio_AMR <- fig1e$AMR/fig1e$Total
fig1e$ratio_AFR <- fig1e$AFR/fig1e$Total
fig1e$ratio_EUR <- fig1e$EUR/fig1e$Total
fig1e$ratio_SAS <- fig1e$SAS/fig1e$Total

fig1e_reaarange <- fig1e[c(5,2,4,6,3,1),]
fig1e_reaarange 
par(mfrow=c(1,5))
pdf("Figure1e_thousand_race.pdf")
barplot(fig1e_reaarange$ratio_EAS,ylim=c(0,0.20),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'EAS' )
barplot(fig1e_reaarange$ratio_AMR,ylim=c(0,0.20),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'AMR' )
barplot(fig1e_reaarange$ratio_AFR,ylim=c(0,0.20),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'AFR' )
barplot(fig1e_reaarange$ratio_EUR,ylim=c(0,0.20),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'EUR' )
barplot(fig1e_reaarange$ratio_SAS,ylim=c(0,0.20),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'SAS' )
dev.off()

### consider Allele frequency 
# fig1e weighted mean
fig1e <- read.table("figure1e_af_mean.txt",header=T)
fig1e <- t(fig1e)
fig1e<- as.data.frame(fig1e)
colnames(fig1e) <- c("Total","Union")
fig1e
fig1e$ratio_union <- fig1e$Union/fig1e$Total/(2504+1823)
fig1e_reaarange <- fig1e[c(5,2,4,6,3,1),]
fig1e_reaarange$ratio_union

library(plotrix)
par(mfrow=c(1,1))
pdf("Figure1e_af_mean.pdf")
#barplot(fig1e_reaarange$ratio_union,ylim=c(0,12),col=c('red','orange','blue','skyblue','green','yellow'),names.arg=c('TSS_coding','TSS_noncoding','Intergenic','Intragenic_coding','Intragenic_noncoding','Not in CGI'),space=1.2,xlab = 'Types of CpG islands', ylab='Mutation rate',border=c('red','orange','blue','skyblue','green','yellow'),cex.names = 0.8, main = 'Union' )
from = 0.003; to = 0.011
gap.barplot(fig1e_reaarange$ratio_union,gap=c(from,to),col=c('red','orange','blue','skyblue','green','yellow'),xlab = 'Types of CpG islands', ylab='Af-weighted mutation rate',main = 'Union' )
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
#axis(2, at=from)
dev.off()

library(psych)

temp_NM_TSS <- read.table("figure1e_af_variation_NM_TSS.txt",header=F)
temp_NR_TSS <- read.table("figure1e_af_variation_NR_TSS.txt",header=F)
temp_intergenic <- read.table("figure1e_af_variation_intergenic.txt",header=F)
temp_NM_intragenic <- read.table("figure1e_af_variation_NM_intragenic.txt",header=F)
temp_NR_intragenic <- read.table("figure1e_af_variation_NR_intragenic.txt",header=F)
temp_noCGI <- read.table("figure1e_af_variation_noCGI.txt",header=F)
summary_af <- rbind(describe(temp_NM_TSS),describe(temp_NR_TSS),describe(temp_intergenic),describe(temp_NM_intragenic),describe(temp_NR_intragenic),describe(temp_noCGI))
length(temp_NM_TSS[temp_NM_TSS>4000,])
mean(temp_NM_TSS[temp_NM_TSS<4000,])
mean(temp_intergenic[temp_intergenic<4000,])
mean(temp_NM_intragenic[temp_NM_intragenic<4000,])

hist(temp_intergenic)
temp_intergenic)

['noCGI', 'NR_TSS', 'NR_intragenic', 'intergenic', 'NM_TSS', 'NM_intragenic', 'miscellaneous']
[51703250, 116304, 73452, 597266, 1982770, 758164, 611056]






### individual box plot
par(mfrow=c(1,3))
fig1e_individual <- read.table("figure1e_individual_union.txt",header=T)
fig1e_individual <- read.table("figure1e_individual_pcawg.txt",header=T)
fig1e_individual <- read.table("figure1e_individual_thousand.txt",header=T)

fig1e_indi_rearrange <- fig1e_individual[,c(5,2,4,6,3,1)]
fig1e_indi_rearrange
pdf("individual_plot.pdf")
boxplot(fig1e_indi_rearrange[-1,]$NM_TSS/fig1e_indi_rearrange[1,]$NM_TSS,
        fig1e_indi_rearrange[-1,]$NR_TSS/fig1e_indi_rearrange[1,]$NR_TSS,
        fig1e_indi_rearrange[-1,]$intergenic/fig1e_indi_rearrange[1,]$intergenic,
        fig1e_indi_rearrange[-1,]$NM_intragenic/fig1e_indi_rearrange[1,]$NM_intragenic,
        fig1e_indi_rearrange[-1,]$NR_intragenic/fig1e_indi_rearrange[1,]$NR_intragenic,
        fig1e_indi_rearrange[-1,]$noCGI/fig1e_indi_rearrange[1,]$noCGI,
        names= c('NM_TSS', 'NR_TSS', 'intergenic', 'NM_intragenic','NR_intragenic', 'noCGI'),
        col= c('red','orange','blue','skyblue','green','yellow'),
        boxwex=0.5,ylab="individual mutation rate",xlab="CpG island types"
        )
stripchart(list(fig1e_indi_rearrange[-1,]$NM_TSS/fig1e_indi_rearrange[1,]$NM_TSS,
                fig1e_indi_rearrange[-1,]$NR_TSS/fig1e_indi_rearrange[1,]$NR_TSS,
                fig1e_indi_rearrange[-1,]$intergenic/fig1e_indi_rearrange[1,]$intergenic,
                fig1e_indi_rearrange[-1,]$NM_intragenic/fig1e_indi_rearrange[1,]$NM_intragenic,
                fig1e_indi_rearrange[-1,]$NR_intragenic/fig1e_indi_rearrange[1,]$NR_intragenic,
                fig1e_indi_rearrange[-1,]$noCGI/fig1e_indi_rearrange[1,]$noCGI
                ), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = c('red','orange','blue','skyblue','green','yellow'))
col = c(alpha('black',0.2),alpha('black',0.2),alpha('black',0.2),alpha('black',0.2),alpha('black',0.2),alpha('black',0.2))

dev.off()



summary(fig1e_indi_rearrange[-1,]$intergenic/fig1e_indi_rearrange[1,]$intergenic)
summary(fig1e_indi_rearrange[-1,]$NM_intragenic/fig1e_indi_rearrange[1,]$NM_intragenic)
summary(fig1e_indi_rearrange[-1,]$NR_intragenic/fig1e_indi_rearrange[1,]$NR_intragenic)


##figure 1f 
fig1f <- read.table("figure1f_af_group.txt",header=T)
fig1f <- t(fig1f)
fig1f<- as.data.frame(fig1f)
colnames(fig1f) <- c("Total","zero","0.06","0.12","0.18","0.24","Other")
fig1f[,c(1:7)]
fig1f$ratio_00 <- fig1f$zero/fig1f$Total
fig1f$ratio_06<- fig1f$`0.06`/fig1f$Total
fig1f$ratio_12 <- fig1f$`0.12`/fig1f$Total
fig1f$ratio_18 <- fig1f$`0.18`/fig1f$Total
fig1f$ratio_24 <- fig1f$`0.24`/fig1f$Total
fig1f$ratio_other <- fig1f$Other/fig1f$Total

fig1f_reaarange <- fig1f[c(5,2,4,6,3,1),]
fig1f_reaarange[,c(1:7)]
fig1f_reaarange_2 <- as.matrix(fig1f_reaarange[,c(8,9,10,11,12,13)])
fig1f_reaarange_2
par(mfrow=c(1,1))
pdf("Figure1f.pdf")
barplot(fig1f_reaarange_2,beside=T,log='y',ylim=c(0.0001,1),col=c('red','orange','blue','skyblue','green','yellow'),xlab='Allele frequency',ylab='Proportion of CpG')
legend(25,0.3,legend=c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI','Not in CGI'),fill=c('red','orange','blue','skyblue','green','yellow'))
dev.off()

###Figure 1b
fib1b <- read.table("CGI_classification_revised.txt",header=F)
fib1b
fig1b_NM_TSS <- fib1b[fib1b$V14=='NM_TSS',]
fig1b_NR_TSS <- fib1b[fib1b$V14=='NR_TSS',]
fig1b_intergenic <- fib1b[fib1b$V14=='intergenic',]
fig1b_NM_intragenic <- fib1b[fib1b$V14=='NM_intragenic',]
fig1b_NR_intragenic <- fib1b[fib1b$V14=='NR_intragenic',]
#fig1b_others <- fib1b[fib1b$V14=='intergenic',]
pdf("fig1b.pdf")
boxplot(fig1b_NM_TSS$V8,fig1b_NR_TSS$V8,fig1b_intergenic$V8,fig1b_NM_intragenic$V8,fig1b_NR_intragenic$V8,outline=F,pars=list(boxwex=0.4,staplewex=0.5,outwex=0.5),names=c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI'),xlab="Type of CpG islands",ylab = 'Size of CpG islands',col=c('red','orange','blue','skyblue','green'))
dev.off()
fig1b_others <- fib1b[-c(fib1b$V14=='intergenic',fib1b$V14=='NR_intragenic',fib1b$V14=='NM_intragenic',fib1b$V14=='NR_TSS',fib1b$V14=='NM_TSS'),]

fig1b_others <- fib1b[fib1b$V14!='intergenic',] %>% .[.$V14!='NM_TSS',] %>% .[.$V14!='NR_TSS',] %>% .[.$V14!='NM_intragenic',] %>% .[.$V14!='NR_intragenic',]
fig1b_others
summary(fig1b_NM_TSS$V8)
summary(fig1b_NR_TSS$V8)
summary(fig1b_intergenic$V8)
summary(fig1b_NM_intragenic$V8)
summary(fig1b_NR_intragenic$V8)
summary(fig1b_others$V8)
summary(fib1b$V8)

length(fig1b_NM_TSS$V8)
length(fig1b_NR_TSS$V8)
length(fig1b_intergenic$V8)
length(fig1b_NM_intragenic$V8)
length(fig1b_NR_intragenic$V8)
length(fig1b_others$V8)
nrow(fib1b)

##Fig 1c
fig1c<-matrix(c(0.48,0.03,0.15,0.18,0.02,0.15),ncol=1)
pdf("fig1c.pdf")
barplot(fig1c,col=c('red','orange','blue','skyblue','green','darkseagreen'))
legend(0.7,0.5,legend=rev(c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI','Miscellaneous')),fill=rev(c('red','orange','blue','skyblue','green','darkseagreen')))
dev.off()

###Figure 2A
fig2<-read.table("distance_by_group_Y.txt",header=F)
colnames(fig2)[1] = 'distance'
fig2
fig2_noCpG<-read.table("../CpG/fig2a_noCpG_mod.txt",header=T)
fig2_noCpG
fig2_new<-as.tibble(fig2) %>% left_join(as.tibble(fig2_noCpG),by = "distance")
tail(fig2_new)
fig2_new

#plot(fig2_new$distance,fig2_new$V3/fig2_new$V2,cex=0.1, xlab='Distance from the nearest CpG island border (base)',ylab='Normalized incidence of C>T mutation',xlim=c(-5500,5500),ylim=c(0.0,0.35),axes=F)
#axis(side=1,at=c(-6000,-5500,-2500,-500,500,2500,5500,6000)) + axis(side=2,at=c(-0.1,0.0,0.1,0.2,0.3,0.4))
#points(fig2_new$distance,fig2_new$incidence,cex=0.1,col='gray')
#dev.off()
#apply(matrix(fig2_new$distance,10),2,mean)
pdf("fig2a.pdf")
plot(apply(matrix(fig2_new$distance,2),2,mean),apply(matrix(fig2_new$V3,2),2,mean)/apply(matrix(fig2_new$V2,2),2,mean),cex=0.1, xlab='Distance from the nearest CpG island border (base)',ylab='Normalized incidence of C>T mutation',xlim=c(-5500,5500),ylim=c(0.0,0.35),axes=F)
axis(side=1,at=c(-6000,-5500,-2500,-500,500,2500,5500,6000)) + axis(side=2,at=c(-0.1,0.0,0.1,0.2,0.3,0.4))
points(apply(matrix(fig2_new$distance,2),2,mean),apply(matrix(fig2_new$incidence,2),2,mean),cex=0.1,col='gray')
dev.off()


?matrix()
#fig2b
fig2b<-read.table("distance_by_group_Y.txt",header=F)
fig2b
pdf("fig2b.pdf")
plot(apply(matrix(fig2b$V1,20),2,mean),apply(matrix(fig2b$V5,20),2,mean)/apply(matrix(fig2b$V4,20),2,mean),cex=0.5,pch=20, col='red',xlab='Distance from the nearest CpG island border (base)',ylab='Normalized incidence of C>T mutation',xlim=c(-5500,5500),ylim=c(0,0.35),axes=F)+
points(apply(matrix(fig2b$V1,100),2,mean),apply(matrix(fig2b$V7,100),2,mean)/apply(matrix(fig2b$V6,100),2,mean),cex=0.5,pch=20,col='orange')+
points(apply(matrix(fig2b$V1,20),2,mean),apply(matrix(fig2b$V9,20),2,mean)/apply(matrix(fig2b$V8,20),2,mean),cex=0.5,pch=20,col='blue')+
points(apply(matrix(fig2b$V1,20),2,mean),apply(matrix(fig2b$V11,20),2,mean)/apply(matrix(fig2b$V10,20),2,mean),cex=0.5,pch=20,col='skyblue')+
points(apply(matrix(fig2b$V1,100),2,mean),apply(matrix(fig2b$V13,100),2,mean)/apply(matrix(fig2b$V12,100),2,mean),cex=0.5,pch=20,col='green')
axis(side=2,at=c(-0.1,0.0,0.1,0.2,0.3,0.4)) +
  axis(side=1,at=c(-6500,-5500,-2500,-500,500,2500,5500,6500))
legend(1400,0.1,legend=c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI'),pch=20,col=c('red','orange','blue','skyblue','green'),cex=0.8)
dev.off()

#'red','orange','blue','skyblue','green'

#fig2c
fig2c<-read.table("methylation_by_distance_Y.txt",header=F)
fig2c
pdf("fig2c.pdf")
plot(apply(matrix(fig2c$V1,20),2,mean),apply(matrix(fig2c$V5,20),2,mean)/apply(matrix(fig2c$V4,20),2,mean),cex=0.5,pch=20, col='red',xlab='Distance from the nearest CpG island border (base)',ylab='Mean methylation proportion',xlim=c(-5500,5500),ylim=c(0,0.9),axes=F)+
points(apply(matrix(fig2c$V1,100),2,mean),apply(matrix(fig2c$V7,100),2,mean)/apply(matrix(fig2c$V6,100),2,mean),cex=0.5,pch=20,col='orange')+
points(apply(matrix(fig2c$V1,20),2,mean),apply(matrix(fig2c$V9,20),2,mean)/apply(matrix(fig2c$V8,20),2,mean),cex=0.5,pch=20,col='blue')+
points(apply(matrix(fig2c$V1,20),2,mean),apply(matrix(fig2c$V11,20),2,mean)/apply(matrix(fig2c$V10,20),2,mean),cex=0.5,pch=20,col='skyblue')+
points(apply(matrix(fig2c$V1,100),2,mean),apply(matrix(fig2c$V13,100),2,mean)/apply(matrix(fig2c$V12,100),2,mean),cex=0.5,pch=20,col='green')
axis(side=2,at=c(-0.1,0.0,0.3,0.6,.9)) +
axis(side=1,at=c(-6500,-5500,-2500,-500,500,2500,5500,6500))
#legend(1400,0.1,legend=c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI'),pch=20,col=c('red','orange','blue','skyblue','green'),cex=0.8)
dev.off()
  
#fig2d
install.packages("vioplot")
library(vioplot)
fig2d<-read.table('methylation_by_CGI.txt',header=T,sep='\t')
fig2d
pdf("fig2d.pdf")
plot(0,0,type='n',xlim=c(0.5,5.5),ylim=c(0,1),axes=F,xlab='CpG island type',ylab='Mean metylation Proportion')
vioplot(fig2d$methylation[fig2d$input_class=='NM_TSS'],at=1,add=T,col='red',border='red',rectCol='red',wex=0.5)
vioplot(fig2d$methylation[fig2d$input_class=='NR_TSS'],at=2,add=T,col='orange',border='orange',rectCol='orange',wex=0.5)
vioplot(fig2d$methylation[fig2d$input_class=='intergenic'],at=3,add=T,col='blue',border='blue',rectCol='blue',wex=0.5)
vioplot(fig2d$methylation[fig2d$input_class=='NM_intragenic'],at=4,add=T,col='skyblue',border='skyblue',rectCol='skyblue',wex=0.5)
vioplot(fig2d$methylation[fig2d$input_class=='NR_intragenic'],at=5,add=T,col='green',border='green',rectCol='green',wex=0.5)
axis(side=1,at=1:5,labels=c('TSS_coding CGI','TSS_noncoding_CGI','Intergenic CGI','Intragenic_coding CGI','Intragenic_noncoding CGI'))+
  axis(side=2,at=c(0,0.33,0.67,1.0))
dev.off()

nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NM_TSS',]) / nrow(fig2d[fig2d$input_class=='NM_TSS',])
nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NR_TSS',]) / nrow(fig2d[fig2d$input_class=='NR_TSS',])
nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='intergenic',]) / nrow(fig2d[fig2d$input_class=='intergenic',])
nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NM_intragenic',]) / nrow(fig2d[fig2d$input_class=='NM_intragenic',])
nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NR_intragenic',]) / nrow(fig2d[fig2d$input_class=='NR_intragenic',])
pdf("fig2d_pie.pdf")
par(mfrow=c(1,5))
pie(c(nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NM_TSS',]),nrow(fig2d[fig2d$methylation<0.67 & fig2d$input_class=='NM_TSS',])),col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,angle=0,radius = 1)+
pie(c(nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NR_TSS',]),nrow(fig2d[fig2d$methylation<0.67 & fig2d$input_class=='NR_TSS',])),col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,angle=0,radius = 1)+
pie(c(nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='intergenic',]),nrow(fig2d[fig2d$methylation<0.67 & fig2d$input_class=='intergenic',])),col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,angle=0,radius = 1)+
pie(c(nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NM_intragenic',]),nrow(fig2d[fig2d$methylation<0.67 & fig2d$input_class=='NM_intragenic',])),col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,angle=0,radius = 1)+
pie(c(nrow(fig2d[fig2d$methylation>=0.67 & fig2d$input_class=='NR_intragenic',]),nrow(fig2d[fig2d$methylation<0.67 & fig2d$input_class=='NR_intragenic',])),col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,angle=0,radius = 1)
dev.off()

fig2d
#fig2e
par(mfrow=c(1,1))
pdf("fig2e.pdf")
boxplot(fig2d$size[fig2d$methylation>=0.67 & fig2d$input_class=='NM_intragenic'],fig2d$size[fig2d$methylation<0.67 & fig2d$input_class=='NM_intragenic'],outline=F,
        pars=list(boxwex=0.3,staplewex=0.5,outwex=0.5),names=c('Group A','Group B'),
        ylab = 'Size of CpG islands',col=c('cyan1','cyan4'),yaxt='n',ylim=c(200,1800))
axis(side=2,at=c(200,600,1000,1400,1800))
dev.off()
#fig2f
sum(fig2d$mutated[fig2d$methylation>=0.67 & fig2d$input_class=='NM_intragenic']) / sum(fig2d$total[fig2d$methylation>=0.67 & fig2d$input_class=='NM_intragenic'])
sum(fig2d$mutated[fig2d$methylation<0.67 & fig2d$input_class=='NM_intragenic'])/sum(fig2d$total[fig2d$methylation<0.67 & fig2d$input_class=='NM_intragenic'])
sum(fig2d$mutated[fig2d$input_class=='NM_intragenic'])/sum(fig2d$total[fig2d$input_class=='NM_intragenic'])
pdf("fig2f.pdf")
barplot(c(0.2630227743,0.05998946654),beside = T,ylim=c(0,0.3),col=c('cyan1','cyan4'),border=c('cyan1','cyan4'),space=2,names.arg = c('Group A', 'Group B'),family='serif',ylab='Mutation rate',xlab='Intragenic_coding CGI')
abline(h=0.276)
abline(h=0.1116695827)
dev.off()

#Fig3
fig3a<-read.table("fig3a_pie.txt",header=T) # methyl>=0.67
fig3a_all<-read.table("fig3a_pie_nomethyl.txt",header=T)
fig3a<-as.matrix(fig3a)
fig3a_all<-as.matrix(fig3a_all)
fig3a[,1]/(fig3a[,1]+fig3a[,2])

pdf("fig3a.pdf")
par(mfrow=c(3,3))
for(i in 1:9) {
  pie(fig3a[i,c(1,2)],col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,radius = 1)
  par(new=T)
  pie(fig3a[i,c(1,2)],col=c('white','white'),labels = NA,border=NA,clockwise = T,radius = 0.8)
  par(new=T)
  pie(fig3a[i,c(3,4)],col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,radius = 0.6)
}
dev.off()

pdf("figs1.pdf")
par(mfrow=c(3,3))
for(i in 1:9) {
  pie(fig3a_all[i,c(1,2)],col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,radius = 1)
  par(new=T)
  pie(fig3a_all[i,c(1,2)],col=c('white','white'),labels = NA,border=NA,clockwise = T,radius = 0.8)
  par(new=T)
  pie(fig3a_all[i,c(3,4)],col=c('violetred1','turquoise3'),labels = NA,border=NA,clockwise = T,radius = 0.6)
}
dev.off()


matrix(fig3a[,1]/(fig3a[,1]+fig3a[,2]),nrow=3,byrow = T)
matrix(fig3a[,3]/(fig3a[,3]+fig3a[,4]),nrow=3,byrow = T)

matrix(fig3a_all[,1]/(fig3a_all[,1]+fig3a_all[,2]),nrow=3,byrow = T)
matrix(fig3a_all[,3]/(fig3a_all[,3]+fig3a_all[,4]),nrow=3,byrow = T)
(2504+1823)

fig3a$total <- fig3a$yesWGS+fig3a$noWGS
fig3a_all$total <- fig3a_all$yesWGS+fig3a_all$noWGS
fig3a
fig3a_all

#Fig3B
fig3b<-read.table("fig3b.txt",header=F) # methyl>=0.67
fig3b
fig3b$total <- apply(fig3b,1,sum)
as.matrix(fig3b[,c(1:6)]/fig3b$total)
pdf("fig3b.pdf")
par(mfrow=c(1,1))
a<-t(as.matrix(fig3b[,c(1:6)]/fig3b$total))
barplot(a[,c(9:1)],beside = F, horiz = T,axes = F,col = rainbow(6))
dev.off()
