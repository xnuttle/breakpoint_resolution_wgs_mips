## makeRunMipDesignShellScript2.py

import sys

regionfile = sys.argv[1]
insertsize = sys.argv[2]
featureflank = sys.argv[3]
mip_overlap = sys.argv[4]
log_prefix = sys.argv[5]
path_to_scripts = sys.argv[6]
genome_dir = sys.argv[7]
fasta_list = sys.argv[8]
check_for_snps=sys.argv[9] #BJO
snp_file=sys.argv[10] #BJO
captureflank=sys.argv[11] #BJO

fnames = []

ear = range(16,21)

#print "rm -rf tmp_dir"

for extarm in ear:
    ligarm = 40-extarm
    logfile = "%s_%d_%d.design_all_mips_joe2" % ( log_prefix, extarm, ligarm )
    #BJO
    cmd = "%s/design_all_mips_joe2.pl -regions_to_scan %s -mip_scan_size %s -extension_arm_length \
%d -ligation_arm_length %d -feature_flank %s -check_for_snps %s -snp_file %s -genome_dir %s > %s.so 2> %s.se &" % \
    ( path_to_scripts,  regionfile, insertsize, extarm, ligarm, featureflank, check_for_snps, snp_file, genome_dir, logfile, logfile )#BKM fixed path
    print "echo \"Running: %s\"" % ( cmd )
    print cmd
    #print "echo \"Finished: %s\"" % ( cmd )
    fnames.append( "%s.ext%d.lig%d.scan%s.all_mips" % ( regionfile, extarm, ligarm, insertsize ) )

print "wait"

for extarm in ear:
    ligarm = 40-extarm
    allmipsfile = "%s.ext%d.lig%d.scan%s.all_mips" % (regionfile, extarm, ligarm, insertsize)
    logfile = "%s_%d_%d.copy_counts" % ( log_prefix, extarm, ligarm )
    cmd = "%s/genome_compare %s %s %s.copy_counts > %s.so 2> %s.se &" % ( path_to_scripts, fasta_list, allmipsfile, allmipsfile, logfile, logfile )
    print "echo \"Running: %s\"" % ( cmd )
    print cmd

print "wait"

fnames = []

for extarm in ear:
    ligarm = 40-extarm
    fnames.append( "%s.ext%d.lig%d.scan%s.all_mips.copy_counts.ranked" % ( regionfile, extarm, ligarm, insertsize ) )
    copycountsfile = "%s.ext%d.lig%d.scan%s.all_mips.copy_counts" % (regionfile, extarm, ligarm, insertsize)
    logfile = "%s_%d_%d.apply_scoring_matrix_in_place4" % ( log_prefix, extarm, ligarm )
    cmd = "python %s/apply_scoring_matrix_in_place4.py %s/scoring_matrix_v4.txt.csv %s %s.ranked %d %d %s > %s.so 2> %s.se" % ( path_to_scripts, path_to_scripts, copycountsfile, copycountsfile, ligarm, extarm, mip_overlap, logfile, logfile )
    print "echo \"Running: %s\"" % ( cmd )
    print cmd
    #print "echo \"Finished: %s\"" % ( cmd )

print "wait"

mips_file_list="%s.scan%s.all_mips.copy_counts.ranked.ranked_list.txt"%(regionfile, insertsize)
f=open( mips_file_list,'w' )
f.write( "\n".join( fnames ) )
f.close()

logfile = "%s.pick_mips4" % ( log_prefix )
#BJO
cmd = "%s/pick_mips4_double_up_on_small.pl -regions_to_scan %s -mips_file_list %s -mip_scan_size %s -mip_overlap %s \
-feature_flank %s -bad_junctions_file %s/bad_junctions.txt > %s.so 2> %s.se" % ( path_to_scripts, regionfile, mips_file_list, insertsize, mip_overlap, captureflank, path_to_scripts, logfile, logfile )
print "echo \"Running: %s\"" % ( cmd )
print cmd
print "echo \"Finished: %s\"" % ( cmd )
