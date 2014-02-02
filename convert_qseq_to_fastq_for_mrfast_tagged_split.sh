if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 <(int)molecular_tag_length>"
    exit 1
fi

if [[ "$1" -gt "4" ]]
then
	J=1
else
  echo "Usage: $0 <(int)molecular_tag_length>"
	echo "Molecular tag length must be at least 5."
  exit 1
fi

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

mkdir ../mrfast_mapping_input/
mkdir ../mrfast_mapping_output/

J=1
for i in $(ls|grep 's_._1_');
  do LANE=$(echo $i|cut -f2 -d_) R1=$i; Rind=$(echo $i|sed "s/[^s]_1_/${LANE}_2_/g"); R2=$(echo $i|sed "s/[^s]_1_/${LANE}_3_/g");
  echo "$BRKPT_SOFTWARE/qseq_to_fastq_for_mrfast_tagged_split $R1 $Rind $R2 76 $1 ${J} 500000 $PARENT_DIR/$BARCODE_KEY;
		mv *_FS${J}_* $PARENT_DIR/mrfast_mapping_input/;"|tr '\n' ' '
	echo
  	J=$(echo "$J+1"|bc)
done > qseq_to_fastq_commands.sh
echo "rm -f qseq_to_fastq_commands.sh" >> qseq_to_fastq_commands.sh

#fastq file naming scheme: <sample name>_FS<file set number>_<part number (corresponding to the 1st, or 2nd, or so on part of the original file set)>.fastq.gz
