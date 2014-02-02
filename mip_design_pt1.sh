#Xander Nuttle
#mip_design_pt1.sh
#Call: bash mip_design_pt1.sh

$BRKPT_SOFTWARE/get_brkpt_SUN_coords prox_aligned.fasta dist_aligned.fasta
$BRKPT_SOFTWARE/brkpt_mipmaker prox_aligned.fasta dist_aligned.fasta
$BRKPT_SOFTWARE/make_SNP_table prox_aligned.fasta dist_aligned.fasta
mkdir genome_for_mip_design
cd genome_for_mip_design
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr*'
md5sum * > md5sums.txt
sort -k2,2 md5sums.txt > dummy1
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/md5*'
sort -k2,2 md5sum.txt > dummy2
diff dummy1 dummy2 > /dev/null
if [ $? -ne 0 ]
	then cd ..
	rm -rf genome_for_mip_design
	echo "DOWNLOAD ERROR! Attempt the download again by rerunning this script."
	exit 1
fi
rm -f dummy*
chmod u+w md5sum.txt
rm -f md5sum*
gunzip *gz
sed 's/prox/chrprox_shared/g' ../prox.fasta > chrprox_shared.fa
for i in `ls *.fa`; do readlink -f $i; done > ../chromosomes.txt
cd ..
python $BRKPT_SOFTWARE/makeRunMipDesignShellScript2.py prox_regions_for_design.bed 112 0 1 brkpt_mips $BRKPT_SOFTWARE ./genome_for_mip_design chromosomes.txt y prox.snptable 0 > makemips.sh
bash makemips.sh &
wait
GENE="prox"
echo "mv ${GENE}_regions_for_design.bed.picked_mip_probe_arms.${GENE}_regions_for_design.bed.scan112.all_mips.copy_counts.ranked.ranked_list.txt ${GENE}.mipdesign"|bash
echo "grep '5\.0' $GENE.mipdesign > ${GENE}_score5.mipdesign"|bash
echo "cut -f10 ${GENE}_score5.mipdesign|sed 's/A/a/g'|sed 's/C/c/g'|sed 's/G/g/g'|sed 's/T/t/g' > ${GENE}_score5.ligarms"|bash
echo "cut -f6 ${GENE}_score5.mipdesign|sed 's/A/a/g'|sed 's/C/c/g'|sed 's/G/g/g'|sed 's/T/t/g' > ${GENE}_score5.extarms"|bash
echo "for i in \`cat ${GENE}_score5.extarms\`; do echo \"CTTCAGCTTCCCGATATCCGACGGTAGTGTNNNNNNNN\" >> ${GENE}_score5.linkers; done"|bash
echo "paste ${GENE}_score5.ligarms ${GENE}_score5.linkers ${GENE}_score5.extarms -d \"\" > ${GENE}_score5.mipoligos"|bash
echo "paste ${GENE}_score5.mipdesign ${GENE}_score5.mipoligos > dummy"|bash
echo "mv dummy ${GENE}_score5.mipdesign"|bash
echo "grep '3\.0' ${GENE}.mipdesign > ${GENE}_score3.mipdesign"|bash
echo "cut -f10 ${GENE}_score3.mipdesign|sed 's/A/a/g'|sed 's/C/c/g'|sed 's/G/g/g'|sed 's/T/t/g' > ${GENE}_score3.ligarms"|bash
echo "cut -f6 ${GENE}_score3.mipdesign|sed 's/A/a/g'|sed 's/C/c/g'|sed 's/G/g/g'|sed 's/T/t/g' > ${GENE}_score3.extarms"|bash
echo "for i in \`cat ${GENE}_score3.extarms\`; do echo \"CTTCAGCTTCCCGATATCCGACGGTAGTGTNNNNNNNN\" >> ${GENE}_score3.linkers; done"|bash
echo "paste ${GENE}_score3.ligarms ${GENE}_score3.linkers ${GENE}_score3.extarms -d \"\" > ${GENE}_score3.mipoligos"|bash
echo "paste ${GENE}_score3.mipdesign ${GENE}_score3.mipoligos > dummy"|bash
echo "mv dummy ${GENE}_score3.mipdesign"|bash
cat prox_score5.mipdesign prox_score3.mipdesign > brkpt.mipdesign
