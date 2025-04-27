setwd
library(multipleROC)
data=read.table("roc.txt", header=T, sep="\t", check.names=F, row.names=1)
a1=multipleROC(group~CSTB,data=data)
a2=multipleROC(group~CSTA,data=data)
a3=multipleROC(group~KRT6B,data=data)
a4=multipleROC(group~KRT16,data=data)
plot<-plot_ROC(list(a1,a2,a3,a4),
               show.eta=FALSE,
               show.sens=FALSE 
)
pdf('roc.pdf',width = 8,height = 8)
plot
dev.off()

