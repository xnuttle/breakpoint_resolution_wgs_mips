#Xander Nuttle
#get_genome.sh
#Call: bash get_genome.sh

mkdir reference_genome
cd reference_genome
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr*'
md5sum * > md5sums.txt
sort -k2,2 md5sums.txt > dummy1
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/md5*'
sort -k2,2 md5sum.txt > dummy2
diff dummy1 dummy2 > /dev/null
if [ $? -ne 0 ]
  then cd ..
  rm -rf reference_genome
  echo "DOWNLOAD ERROR! Attempt the download again by rerunning this script."
  exit 1
fi
rm -f dummy*
chmod u+w md5sum.txt
rm -f md5sum*
zcat * >> hg19_unmasked.fasta
cd ..
