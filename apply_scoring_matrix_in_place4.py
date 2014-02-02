import sys

from uniformity_functions import *

if len(sys.argv) != 7:
	print "Usage: python apply_scoring_matrix.py <scoring_matrix> <exome_defs> <out_file> <ligation_length> <extension_length> <bin_size>"
	exit

score_file = sys.argv[1]
exome_defs_file = sys.argv[2]
out_file = sys.argv[3]
ligation_length = int(sys.argv[4])
extension_length = int(sys.argv[5])
bin_size = int(sys.argv[6])

#populate the scores dictionary
fin = open(score_file)
score_raw_data = fin.readlines()
fin.close()
score_raw_data = [entry.strip() for entry in score_raw_data]
scores_dict = {}
for entry in score_raw_data:
	tm_ligation = float(entry.split('\t')[0])
	tm_extension = float(entry.split('\t')[1])
	score = float(entry.split('\t')[2])
	key = (tm_ligation, tm_extension)
	scores_dict[key] = score

#initialize the in and out files
fin = open(exome_defs_file)
fout = open(out_file,'w')
for i in range(6):
	fout.write(fin.readline())
fin.readline()
new_header = ">mip_count	[rank_score]	chr	ext_probe_start ext_probe_stop	ext_probe_sequence	ext_probe_copy_count	lig_probe_start lig_probe_stop	lig_probe_sequence	lig_probe_copy_count	mip_scan_start_position	mip_scan_stop_position	scan_target_sequence	feature_start_position	feature_stop_position	feature_mip_count	probe_strand	notes"
fout.write(new_header + '\n')
for line in fin:
	values = line.strip().split('\t')
	if len(values) > 4:
		mip_index = values[0]
		chr = values[1]
		ext_probe_start = values[2]
		ext_probe_stop = values[3]
		extension_sequence = values[4][0:extension_length]
		ext_probe_copy_count = values[5]
		lig_probe_start = values[6]
		lig_probe_stop = values[7]
		ligation_sequence = values[8][0:ligation_length]
		lig_probe_copy_count = values[9]
		mip_scan_start_position = values[10]
		mip_scan_stop_position = values[11]
		scan_target_sequence = values[12]
		feature_start_position = values[13]
		feature_stop_position = values[14]
		feature_mip_count = values[15]
		probe_strand = values[16]
		if (len(values)>17):
			notes = values[17]
		else:
			notes = ""
		# remove indels BJO
		
		if 'L' in notes or '-' in notes or 'D' in notes or '(' in notes or '(/L' in notes or 'l/e' in notes:
			pass
		else:
			#print notes
			if (extension_sequence.count('n')==0 and ligation_sequence.count('n')==0 and extension_sequence.count('i')==0 and ligation_sequence.count('i')==0 and extension_sequence.count('-')==0 and ligation_sequence.count('-')==0 and extension_sequence.count('N') == 0 and ligation_sequence.count('N')==0):
			#valid MIP to calculate Tm
				extension_tm = calc_tm_sugi([extension_sequence], 0.72*pow(10,-12), 0.1)
				extension_tm = extension_tm[0]
				ligation_tm = calc_tm_sugi([ligation_sequence], 0.72*pow(10,-12), 0.1)
				ligation_tm = ligation_tm[0]
				tm_lig_bin = round(ligation_tm/float(bin_size)) * bin_size
				tm_ext_bin = round(extension_tm/float(bin_size)) * bin_size
				key = (tm_lig_bin, tm_ext_bin)
				if scores_dict.has_key(key):
					fout.write(mip_index + '\t' + str(scores_dict[key]) + '\t' + chr + '\t' + ext_probe_start + '\t' + ext_probe_stop + '\t' + extension_sequence + '\t' + ext_probe_copy_count + '\t' + lig_probe_start + '\t' + lig_probe_stop + '\t' + ligation_sequence + '\t' + lig_probe_copy_count + '\t' + mip_scan_start_position + '\t' + mip_scan_stop_position + '\t' + scan_target_sequence + '\t' + feature_start_position + '\t' + feature_stop_position + '\t' + feature_mip_count + '\t' + probe_strand + '\t' + notes +  '\n')
				else:
					fout.write(mip_index + '\t' + "-1" + '\t' + chr + '\t' + ext_probe_start + '\t' + ext_probe_stop + '\t' + extension_sequence + '\t' + ext_probe_copy_count + '\t' + lig_probe_start + '\t' + lig_probe_stop + '\t' + ligation_sequence + '\t' + lig_probe_copy_count + '\t' + mip_scan_start_position + '\t' + mip_scan_stop_position + '\t' + scan_target_sequence + '\t' + feature_start_position + '\t' + feature_stop_position + '\t' + feature_mip_count + '\t' + probe_strand + '\t' + notes +  '\n')
fin.close()
fout.close()
