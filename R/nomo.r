setwd
library(rms)
library(pROC)
library(RColorBrewer)
group <- read.csv("group.csv",header=T,row.names = 1,check.names = F,stringsAsFactors = F)
exp <- read.csv("biomarker.csv",header=T,row.names = 1,check.names = F,stringsAsFactors = F)

exp <-as.data.frame(exp)
f_lrm <-lrm(group ~ ., data=exp, x=TRUE, y=TRUE, maxit=1000)

ddist <- datadist(exp)
options(datadist='ddist')

nomogram <- nomogram(f_lrm,fun=function(x)1/(1+exp(-x)), 
                     fun.at = c(0.1,0.5,0.999),
                     funlabel = "Probability of Case", 
                     lp=F,  
                     conf.int = F, 
                     abbrev = F
)

pdf(file="01_nomogram.pdf" ,width=12,height=8) 

plot(nomogram, lplabel = "Linear Predictor", xfrac = 0.1, varname.label = TRUE, varname.label.sep = "=",
     ia.space = 0.2, tck = NA, tcl = -0.2, lmgp = 0.3, points.label = "Points", total.points.label = "Total Points",
     total.sep.page = FALSE, cap.labels = FALSE, cex.var = 1, cex.axis = 0.6, lwd = 5,
     label.every = 1, col.grid = gray(c(0.8, 0.95)))

dev.off()

exp_logic <- data.frame(group = group, exp)
f_lrm <-lrm(group ~ ., data=exp_logic, x=TRUE, y=TRUE, maxit=1000)

set.seed(123)
cal <- calibrate(f_lrm)
pdf(file="03_calibrate.pdf",width=5,height=5)
par(mgp = c(2.5, 1, 0), mar = c(5.5, 4, 1, 1))
plot(cal)
dev.off()

prob <-  predict(f_lrm, exp_logic[, -1])
pred_LR <- data.frame(prob, group = group, stringsAsFactors = F)

pdf(file="02_roc.pdf" ,width=4.5,height=4.5) 
LR.roc <- plot.roc(pred_LR$group,pred_LR$prob,ylim=c(0,1),xlim=c(1,0),
                   smooth=F,
                   ci=TRUE,
                   main="",
                   lwd=2, 
                   legacy.axes=T,
                   print.auc.col="red",
                   print.auc.x=0.7,
                   print.auc.y=0.3,
                   print.auc = T)
dev.off()


library(rmda)

exp_logic$group <- ifelse(exp_logic$group == "Control", 0, 1)
Nomogram <- decision_curve(group ~ CSTB + CSTA + KRT16 , data = exp_logic
                           
                           ,study.design = 'cohort')
CSTB<- decision_curve(group ~ CSTB, data = exp_logic
                     
                     ,study.design = 'cohort')
CSTA <- decision_curve(group ~ CSTA , data = exp_logic
                        
                        ,study.design = 'cohort')
KRT16 <- decision_curve(group ~ KRT16, data = exp_logic
                       
                       ,study.design = 'cohort')


list <- list(Nomogram, CSTB, CSTA, KRT16)

pdf("04_DCA.pdf",width=5,height=5)
par(mgp = c(2.5, 1, 0), mar = c(4, 4.5, 1, 1))
plot_decision_curve(list,
                    curve.names= c("Nomogram", "CSTB", "CSTA", "KRT16"),
                    cost.benefit.axis =FALSE,
                    col= brewer.pal(6, "Set3"),
                    confidence.intervals=FALSE,
                    standardize = FALSE,
                    legend.position="bottomleft")
dev.off()
