#Xander Nuttle
#mip_analysis_pt1.sh
#Call: bash mip_analysis_pt1.sh

bash $BRKPT_SOFTWARE/mrzip.sh
bash zip_commands.sh
bash $BRKPT_SOFTWARE/convert_qseq_to_fastq_for_mrfast_tagged_split.sh 8
bash qseq_to_fastq_commands.sh
