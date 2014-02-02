#Xander Nuttle
#wgs_analysis_pt2.sh
#Call: bash wgs_analysis_pt2.sh

bash $BRKPT_SOFTWARE/mrzip.sh
bash zip_commands.sh
bash $BRKPT_SOFTWARE/convert_qseq_to_fastq_for_mrfast_split.sh
