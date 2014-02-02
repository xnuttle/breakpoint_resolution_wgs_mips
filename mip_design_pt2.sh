#Xander Nuttle
#mip_design_pt2.sh
#Call: bash mip_design_pt2.sh

$BRKPT_SOFTWARE/pickmips brkpt.filtered.mipdesign brkpt.suns 2
cut -f3 brkpt.filtered.mippicks > dummy1
cut -f4-5 brkpt.filtered.mippicks > dummy2
cut -f8-9 brkpt.filtered.mippicks > dummy3
paste dummy1 dummy2 dummy3 > brkpt_mip.armlocs
sed 's/prox_shared/prox/g' brkpt_mip.armlocs > dummy
mv dummy brkpt_mip.armlocs
rm -f dummy*
cut -f1-4 brkpt.suns > brkpt.suns.reduced
echo "prox" > dummy1
echo "0" > dummy2
echo "0" > dummy3
paste dummy1 dummy2 dummy3 > dummy.exons
$BRKPT_SOFTWARE/detail_mip_targets brkpt_mip.armlocs 1 2 prox.fasta prox.fasta dist.fasta brkpt.suns.reduced dummy.exons brkpt
rm -f dummy*
echo "MIP design complete! Final designed MIPs are in 'brkpt.filtered.mippicks'."
echo "The file 'brkpt.miptargets' has also been generated and will be needed during data analysis."
