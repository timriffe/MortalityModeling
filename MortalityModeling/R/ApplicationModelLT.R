
# Author: tim
###############################################################################
setwd("/home/tim/workspace/MortalityModeling")
#library(devtools)
#install_github("mpascariu/MortalityEstimate",force=TRUE)
library(MortalityEstimate)
# colors. 
library(RColorBrewer)

# subset mortality from which to derive standard
HMD719m <- HMD719[HMD719$sex == "male", ]

# lower age bounds
x <- c(0,1, seq(5, 110, by = 5))

# derive standard (takes ~20 seconds)
W <- MortalityEstimate::wilmoth(x, LT = HMD719m)

# select different combos of child and adult mort
rat    <- c(1.5,2.5,5)
q0_5   <- c(.01,.05,.1)
q15_45 <- outer(q0_5,rat,"*")

# some colors to vary
reds    <- brewer.pal(5,"Reds")[3:5]
blues   <- brewer.pal(5,"Blues")[3:5]
purples <- brewer.pal(5,"Purples")[3:5]
CC      <- rbind(blues,purples,reds)

LX <- matrix(ncol = 9, nrow = length(x))
for (i in 1:length(q0_5)){
	for (j in 1:ncol(q15_45)){
		LX[,(i*3)-3+j] <- wilmothLT(W, q0_5 = q0_5[i], q15_45 = q15_45[i,j], radix = 1)$lt$lx
	}
}

# need to manually add LX <-  to top of this output
dput(LX,file = "Data/LX.R")
#
#matplot(x, LX, lwd = 2, col = t(CC), lty = 1:3,type='l',
#		xlab="Age",ylab="Survival probability", las = 1)
#text(3,LX[3, c(1,4,7)]+.015, c("0.01","0.05","0.10"), pos = 4)
#text(0,.7,"child mortality\n(ages 0-5)",pos=4)
#arrows(4,.76,8,.88, length = .08)
#
#text(28,LX[8, c(7:9)]+ c(.01,-.03,-.06), c("1.5x","2.5x","5.0x") , pos = 4)
#text(20,.5,"adult mortality\n(ages 15-65)\nrelative to child",pos=4)
#arrows(30,.6,32,.7, length = .08)
#
#locator(2)
#
#LX[3, c(1,4,7)]
## save out simple plot, then mark it up in Inkscape,
## then save out as 1200 dpi eps file
##pdf("Figs/WilmothVariants.pdf")
##plot(NULL,type='l',xlim=c(0,110),ylim=c(0,1),xlab="Age",ylab="Survivorship", las = 1)
##for (i in 1:length(q0_5)){
##	for (j in 1:ncol(q15_45)){
##		lines(x, wilmothLT(W, q0_5 = q0_5[i], q15_45 = q15_45[i,j], radix = 1)$lt$lx,
##				col = CC[i,j], lty=j,lwd=2)
##	}
##}
##dev.off()

# 
# plot(NULL,type='l',xlim=c(0,110),ylim=c(1e-5,1),xlab="Age",ylab="log(mx)",log='y')
# for (i in 1:length(q0_5)){
# 	for (j in 1:ncol(q15_45)){
# 		lines(x, wilmothLT(W, q0_5 = q0_5[i], q15_45 = q15_45[i,j])$lt$mx,col="#00000050")
# 	}
# }

