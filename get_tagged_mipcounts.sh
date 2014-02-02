if [[ "$#" -ne "2" ]]
then
    echo "Usage: $0 <path_to_miptargets_file/miptargets_file> <(int)molecular_tag_length>"
    exit 1
fi

num_barcodekey_files=$(ls ..|grep barcodekey|wc -l)
if [[ "$num_barcodekey_files" -gt "1" ]]
then
	echo "There are multiple barcodekey files. Remove those that are not applicable for the data here!"
	exit 1
fi

CURRENT_DIR=$(pwd)
PARENT_DIR=$(dirname $CURRENT_DIR)
EXPERIMENT_NAME=$(basename $PARENT_DIR)

barcodefile=$(ls ..|grep barcodekey)
J=1
for i in $(cut -f1 ../$barcodefile);
	do ls|grep $i > files.txt;
	echo "$BRKPT_SOFTWARE/tagged_mrfast_output_to_mipcounts $1 files.txt $2 ${J} >> $EXPERIMENT_NAME.mipcounts"|bash;
	J=$(echo "$J+1"|bc);
done

