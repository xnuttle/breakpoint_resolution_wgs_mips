experiment_name<-paste(base_name,".mipcounts",sep="")
barcode_file<-paste(base_name,".barcodekey",sep="")
key<-read.table(barcode_file,header=F,sep='\t',colClasses="character")
pdfname<-paste(base_name,".pdf",sep="")
pdf(pdfname,width=11,height=7,useDingbats=FALSE)

for(indiv in 1:length(key[,1]))
{
	individual<-key[indiv,1]
	source("mipplot_brkpt.r")
	p<-p+opts(title=individual)
	print(p)
}

dev.off()

