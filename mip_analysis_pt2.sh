#Xander Nuttle
#mip_analysis_pt2.sh
#Call: bash mip_analysis_pt2.sh
#
#run this script from an experiment's "mrfast_mapping_output" directory

rm -rf ../mrfast_mapping_input
rm -f *DIVET*
rm -f *OEA*
rm -f *unmapped*
bash $BRKPT_SOFTWARE/get_tagged_mipcounts.sh ../brkpt.miptargets 8 &
wait
echo "Analysis complete!"
