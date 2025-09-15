pdf("seq3.pdf", width=12, heigh=6)
data<-read.table("seq3.csv", header=T);
maxScale<-1.2*max(c(max(data$Disrupted_Intersect_Region)), max(data$count))
plot(data$count ~ data$Position,type="l",xlab="Position", ylab="Count", ylim=c(0,maxScale))
lines(data$Disrupted_Intersect_Region ~ data$Position,col="red")
lines(data$open ~ data$Position,col="blue")

legend("topleft", legend=c("Count of disruptive SNP in the last 10 nts", "Count of regions significantly disrupted by SNPs", "Opening Energy"),col=c("black","red","blue"),lty=1)

dev.off()
