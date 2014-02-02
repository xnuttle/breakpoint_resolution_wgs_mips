library(ggplot2)
data<-read.table(experiment_name,header=T)
rows<-which(data[,1]==individual)
data2<-data[1:length(rows),]
for(i in 1:length(rows))
{
	data2[i,]<-data[rows[i],]
}
rows<-which(data2[,4]!="E")
data5<-data2[1:length(rows),]
for(i in 1:length(rows))
{
	data5[i,]<-data2[rows[i],]
}
x<-data5[,3]
f1<-x
f2<-x
for(i in 1:length(x))
{
	f1[i]<-(data5[i,5]/(sum(data5[i,5:6],na.rm=TRUE)))
	f2[i]<-(data5[i,6]/(sum(data5[i,5:6],na.rm=TRUE)))
}
freqs<-data.frame(x,f1,f2)
p<-ggplot(freqs,aes(x,f1))
p<-p+geom_point(color="#EE0000")
p<-p+geom_point(data=freqs,aes(x,f2),color="#0000FF")
p<-p+opts(panel.background = theme_rect(fill='darkgray'))
p<-p+coord_cartesian(ylim=c(-0.25,1.1))
p<-p+scale_x_continuous(expression("Alignment Coordinate"))
p<-p+ylab("Paralog-Specific Count Frequency")
p<-p+opts(legend.position="none")
p<-p+opts(title=individual)
p<-p+opts(panel.grid.major=theme_line(colour="lightgray"))
p<-p+opts(panel.grid.minor=theme_line(colour="lightgray"))

