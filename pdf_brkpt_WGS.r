pdfname<-paste(indiv,"_brkpt.pdf",sep="")
pdf(pdfname,width=11,height=7,useDingbats=FALSE)
library(ggplot2)
data<-read.table("brkpt.suns.depth",header=F)
seqs<-data[,1]
aligncoords<-data[,3]
sundepths<-data[,4]
finaldata<-data.frame(seqs,aligncoords,sundepths)
plot<-ggplot(finaldata,aes(aligncoords,sundepths))+geom_bar(color="#000000",size=1,stat="identity")+facet_grid(seqs ~ .)
plot<-plot+opts(panel.background = theme_rect(fill='white'))
plot<-plot+scale_x_continuous(expression("Alignment Coordinate"))
plot<-plot+ylab("Paralog-Specific Read Depth")
plot<-plot+opts(title=indiv)
plot<-plot+opts(legend.position="none")
plot<-plot+opts(panel.grid.major=theme_line(colour="#DDDDDD"))
print(plot)
dev.off()

