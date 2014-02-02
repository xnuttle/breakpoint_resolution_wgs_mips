#Xander Nuttle
#wgs_analysis_pt3.sh
#Call: bash wgs_analysis_pt3.sh

barcodefile=$(ls ..|grep barcodekey)
CURRENT_DIR=$(pwd)
PARENT_DIR=$(dirname $CURRENT_DIR)
BARCODE_KEY=$(basename $barcodefile)

mv *fastq.gz ../fastqs
cd ../fastqs
for i in $(cut -f1 $PARENT_DIR/$BARCODE_KEY); do mkdir $i; mv $i*fastq.gz $i; done
for i in $(cut -f1 $PARENT_DIR/$BARCODE_KEY);
	do cd $i;
	ls|grep fastq > files.txt;
	$BRKPT_SOFTWARE/get_sunk_reads ../../brkpt.sunsunks files.txt;
	rm -f files.txt;
	$BRKPT_SOFTWARE/sundepth ../../brkpt.suns ../../brkpt.sunsunks filtered_reads.fastq;
	cd ..;
done
