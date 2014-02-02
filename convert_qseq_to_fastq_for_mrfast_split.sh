num_barcodekey_files=$(ls ..|grep barcodekey|wc -l)
if [[ "$num_barcodekey_files" -gt "1" ]]
then
  echo "There are multiple barcodekey files. Remove those that are not applicable for the data here!"
  exit 1
fi

barcodefile=$(ls ..|grep barcodekey)
CURRENT_DIR=$(pwd)
PARENT_DIR=$(dirname $CURRENT_DIR)
BARCODE_KEY=$(basename $barcodefile)

mkdir ../fastqs/
ls|grep 's_._1_' > files1
ls|grep 's_._2_' > files2
ls|grep 's_._3_' > files3
paste files1 files2 files3 > qseqs.txt
rm -f files*
$BRKPT_SOFTWARE/qseq_to_fastq_for_mrfast_split qseqs.txt 100 500000 $PARENT_DIR/$BARCODE_KEY
mv *fastq.gz ../fastqs
rm -f qseqs.txt

#fastq file naming scheme: <sample name>_<part number (where each part contains up to 500,000 reads)>.fastq.gz
