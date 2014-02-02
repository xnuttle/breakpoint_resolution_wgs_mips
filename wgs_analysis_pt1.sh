#Xander Nuttle
#wgs_analysis_pt1.sh
#Call: bash wgs_analysis_pt1.sh

$BRKPT_SOFTWARE/get_brkpt_SUN_coords prox_aligned.fasta dist_aligned.fasta
bash $BRKPT_SOFTWARE/get_genome.sh
$BRKPT_SOFTWARE/sunksearch reference_genome/hg19_unmasked.fasta prox.fasta dist.fasta seqs.refcoords brkpt.sunks
if [ $? -ne 0 ]
  then exit 1
fi
$BRKPT_SOFTWARE/sunk_sun_intersect brkpt.suns brkpt.sunks
