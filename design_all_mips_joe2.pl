#!/usr/bin/perl
######################################################
#written by Ian Byrell Stanaway bard@u.washington.edu#
######################################################
use strict;
my %arg;
&parseCommandLine;	
my $regions_to_scan = $arg{-regions_to_scan};
my @regions_to_scan_contents;
my $feature_start_position_flanked;
my $feature_stop_position_flanked;
my $feature_start_position;
my $feature_stop_position;
&get_regions_to_scan;
my $extension_arm_length = $arg{-extension_arm_length};
my $ligation_arm_length = $arg{-ligation_arm_length};
my $mip_scan_size = $arg{-mip_scan_size};
my $feature_flank = 0;
if (exists $arg{-feature_flank})
{
	$feature_flank = $arg{-feature_flank};
}
my $check_for_snps = $arg{-check_for_snps};
my $genome_sequence_dir = $arg{-genome_dir};
my $snp_file = $arg{-snp_file};

my $ext_probe_sequence_rev;
my $lig_probe_sequence_rev;
my $ext_probe_sequence_fwd;
my $lig_probe_sequence_fwd;
my $scan_target_sequence;

#restriction sites
my $ntalwl_seq_5 = "GGATC";
my $ntalwl_seq_3 = "GATCC";
my $nbbsrdi_seq_5 = "GCAATG";
my $nbbsrdi_seq_3 = "CATTGC";

my $bad_resriction_site_design_count = 0;
my $universal_left_end_mip_seq = "AGGACCGGATCAACT";#matches $ntalwl_seq_5 = "GGATC"; once
my $universal_right_end_mip_seq = "CATTGCGTGAACCGA";#matches $nbbsrdi_seq_3 = "CATTGC"; once
my $universal_middle_mip_seq = "CTTCAGCTTCCCGATATCCGACGGTAGTGT";
#example:
#AGGACCGGATCAACTacgcgtgccatctgccacccCTTCAGCTTCCCGATATCCGACGGTAGTGTcaacactcgttttgtgtcccCATTGCGTGAACCGA

my $ext_probe_start_fwd;
my $ext_probe_stop_fwd;
my $lig_probe_start_fwd;
my $lig_probe_stop_fwd;
my $ext_probe_start_rev;
my $ext_probe_stop_rev;
my $lig_probe_start_rev;
my $lig_probe_stop_rev;
my %mip_ext_lig_seq_chr;
my $feature_mip_count = 0;
my $probe_copy_number = 0;
my $mip_count = 0;
my $snp_position;
my $total_snp_load_count = 0;
my $allele_1;
my $allele_2;
my $mip_scan_start_position;
my $mip_scan_stop_position;
my $mip_finished;# y/n
my $design_attempts = 0;
my $design_done = 0;
my $chr_fasta = "NA";
my $chr_sequence = "";
my $chr;
my $query_start;
my $query_stop;
my %chr_snp_positions;
&load_chr_snps;

##############my $genome_sequence_dir = "/net/nickerson/data/ucscHg18/";
#my $genome_sequence_dir = "/droog/ucscHg18/";
#my $genome_sequence_dir = "/nfs/home/bard/ucsc/hg18/";
##############my $command = "mkdir tmp_dir && cp $genome_sequence_dir*.fa ./tmp_dir/";
##############system("$command");
##############print "copied chr fastas to tmp_dir\n";
open (CHRSNPFILE, ">Chr_snps_size.txt") or die "can't create chr_snps_size\n";
my $outfile = "$regions_to_scan" . ".ext$extension_arm_length.lig$ligation_arm_length.scan$mip_scan_size.all_mips";
open (OUTFILE, ">$outfile");
print OUTFILE ">regions_to_scan:$regions_to_scan\n";
print OUTFILE ">extension_arm_length:$extension_arm_length\n";
print OUTFILE ">ligation_arm_length:$ligation_arm_length\n";
print OUTFILE ">mip_scan_size:$mip_scan_size\n";
print OUTFILE ">feature_flank:$feature_flank\n";
print OUTFILE ">max_snp:1\n";
print OUTFILE ">mip_count\tchr\text_probe_start\text_probe_stop\text_probe_sequence\tlig_probe_start\tlig_probe_stop\tlig_probe_sequence\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tfeature_start_position\tfeature_stop_position\tfeature_mip_count\tprobe_strand\tnotes\n";
my $line_count = 0;
my $number_of_features = $#regions_to_scan_contents;
foreach my $lines (@regions_to_scan_contents)
{
	chomp $lines;
	if ($lines =~ m/^>|^#/g)#header
	{}
	else
	{
		$line_count++;
		#print $lines . "\n";
		my @line_contents = split(/\s+/, $lines);
		my $tmp_chr = $chr;#save the old one
		$chr = $line_contents[0];
		$chr =~ s/chr//g;
		$feature_start_position = $line_contents[1] + 1;#make 1 based
		$feature_stop_position = $line_contents[2];
		$feature_start_position_flanked = $line_contents[1] - $feature_flank;
		$feature_stop_position_flanked = $line_contents[2] + $feature_flank;
		
		############## COMMENTED OUT BY JOE ON 110613
		# my $feature_length = $feature_stop_position_flanked - $feature_start_position_flanked + 1;
		# my $tmp_chr_fasta = "./tmp_dir/chr$chr" . ".fa";
		# if ($tmp_chr_fasta ne $chr_fasta)
		# {
			# $chr_fasta = "./tmp_dir/chr$chr" . ".fa";
			# &get_chr_fasta_sequence;
		
		############## ADDED BY JOE ON 110613
		my $tmp_chr_fasta = "$genome_sequence_dir/chr$chr" . ".fa";
		if ($tmp_chr_fasta ne $chr_fasta)
		{
			$chr_fasta = "$genome_sequence_dir/chr$chr" . ".fa";
			&get_chr_fasta_sequence;
			
		# NOT COMMENTED OUT BY JOE:
		#	if ($check_for_snps eq "y")
		#	{
		#		print "loading chr$chr snps\n";
		#		delete ($chr_snp_positions{$tmp_chr});
		#		&load_chr_snps;
		#	}
		}
		$feature_mip_count = 0;
		$mip_scan_stop_position = 0;
		while ($mip_scan_stop_position <= $feature_stop_position_flanked)
		{
			$feature_mip_count++;
			if ($feature_mip_count > 1)
			{
				$mip_scan_start_position++;
				$mip_scan_stop_position++;

				#fwd/+ strand
				$ext_probe_start_fwd++;
				$ext_probe_stop_fwd++;
				$lig_probe_start_fwd++;
				$lig_probe_stop_fwd++;

				#rev/- strand
				$ext_probe_start_rev++;
				$ext_probe_stop_rev++;
				$lig_probe_start_rev++;
				$lig_probe_stop_rev++;
			}
			else
			{
				$mip_scan_start_position = $feature_start_position_flanked;
				$mip_scan_stop_position = $mip_scan_start_position + $mip_scan_size - 1;

				#rev/- strand probes
				$ext_probe_start_rev = $mip_scan_stop_position + 1;
				$ext_probe_stop_rev = $mip_scan_stop_position + $extension_arm_length;
				$lig_probe_start_rev = $mip_scan_start_position - $ligation_arm_length;
				$lig_probe_stop_rev = $mip_scan_start_position - 1;

				#fwd/+ strand probes
				$lig_probe_start_fwd = $mip_scan_stop_position + 1;
				$lig_probe_stop_fwd = $mip_scan_stop_position + $ligation_arm_length;
				$ext_probe_start_fwd = $mip_scan_start_position - $extension_arm_length;
				$ext_probe_stop_fwd = $mip_scan_start_position - 1;
			}
			$scan_target_sequence = uc(substr($chr_sequence,$mip_scan_start_position - 1,$mip_scan_size));		
			&design_mip;
		}
	}
}
close OUTFILE;
print "loaded total $total_snp_load_count snps\nbad_resriction_site_design_count\t$bad_resriction_site_design_count\n";
print CHRSNPFILE "loaded total $total_snp_load_count snps\nbad_resriction_site_design_count\t$bad_resriction_site_design_count\n";
close CHRSNPFILE;

#$command = "rm -r tmp_dir";
#system("$command");

sub design_mip
{
	my $snp_count_fwd = 0;
	for my $int ($ext_probe_start_fwd .. $ext_probe_stop_fwd)
	{
		if (exists $chr_snp_positions{$chr}{$int})
		{
			$snp_count_fwd++;
			$snp_position = $int;
		}
	}
	for my $int ($lig_probe_start_fwd .. $lig_probe_stop_fwd)
	{
		if (exists $chr_snp_positions{$chr}{$int})
		{
			$snp_count_fwd++;
			$snp_position = $int;
		}
	}
	$ext_probe_sequence_fwd = uc(substr($chr_sequence,$ext_probe_start_fwd - 1,$extension_arm_length));
	$lig_probe_sequence_fwd = uc(substr($chr_sequence,$lig_probe_start_fwd - 1,$ligation_arm_length));

	#check forward for restriction sites
	my $ntalwl_seq_5_count = 0;
	my $ntalwl_seq_3_count = 0;
	my $nbbsrdi_seq_5_count = 0;
	my $nbbsrdi_seq_3_count = 0;
	my $scan_target_sequence_with_probe_arms = $ext_probe_sequence_fwd . $scan_target_sequence . $lig_probe_sequence_fwd;
	#print "$scan_target_sequence_with_probe_arms\n";
	#CCTCCCCTGGGCCCCCCACCGGCACCCTCCGCCGCCCCTTCTTGAACACACTCAGGAATGGCTGCTCCAGGTTTTCTCGCTGGTTCTCCAGGTCCAGGATGTCCTAGGAGGAGTAGAGCTCAGGGGAGGGGGCTTTTCCCAGGTTCCTCACA
	my $mip_seq = $universal_left_end_mip_seq . $lig_probe_sequence_fwd . $universal_middle_mip_seq . $ext_probe_sequence_fwd . $universal_right_end_mip_seq;
	#print "$mip_seq\n";
	#AGGACCGGATCAACTCTTTTCCCAGGTTCCTCACACTTCAGCTTCCCGATATCCGACGGTAGTGTCCTCCCCTGGGCCCCCCACCCATTGCGTGAACCGA
	#check the mip for restriction sites
	my $length_mip_seq = length($mip_seq) - 1;#0 based
	for my $int (0 .. ($length_mip_seq - 4))#5 bp window
	{
		my $sub_mip_seq_length_5 = substr($mip_seq, $int, 5);
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_mip_seq - 5))#6 bp window
	{
		my $sub_mip_seq_length_6 = substr($mip_seq, $int, 6);
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}

	#check the target with probe arms
	my $length_target_sequence_with_probe_arms = length($scan_target_sequence_with_probe_arms) - 1;#0 based
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 4))
	{
		my $sub_target_seq_length_5 = substr($length_target_sequence_with_probe_arms, $int, 5);
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 5))
	{
		my $sub_target_seq_length_6 = substr($length_target_sequence_with_probe_arms, $int, 6);
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}

	if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
	{
		if ($snp_count_fwd == 0)
		{
			$mip_count++;
			print OUTFILE "$mip_count\t$chr\t$ext_probe_start_fwd\t$ext_probe_stop_fwd\t$ext_probe_sequence_fwd\t$lig_probe_start_fwd\t$lig_probe_stop_fwd\t$lig_probe_sequence_fwd\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\n";
			print "$line_count/$number_of_features\t$mip_count\n";
		}
		elsif ($snp_count_fwd == 1)
		{#make for allele_1
			#make two mips, one on each strand for each allele of the snp single snp found
			&get_alleles_for_position;#$snp_position
			if ($snp_position >= $ext_probe_start_fwd && $snp_position <= $ext_probe_stop_fwd)
			{
				$mip_count++;
				my $allele_probe_seq_position = $snp_position - $ext_probe_start_fwd;#0 based for array
				my @arr_ext_probe_sequence =  split(//, $ext_probe_sequence_fwd);
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$ext_probe_sequence_fwd = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_sequence_fwd = $ext_probe_sequence_fwd . $bases;
				}
				my $allele_probe_position = $allele_probe_seq_position + 1;
				print OUTFILE "$mip_count\t$chr\t$ext_probe_start_fwd\t$ext_probe_stop_fwd\t$ext_probe_sequence_fwd\t$lig_probe_start_fwd\t$lig_probe_stop_fwd\t$lig_probe_sequence_fwd\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
				print "$line_count/$number_of_features\t$mip_count\n";
			}
			elsif ($snp_position >= $lig_probe_start_fwd && $snp_position <= $lig_probe_stop_fwd)						
			{#must be between the ligerse target arms
				$mip_count++;
				my $allele_probe_seq_position = $snp_position - $lig_probe_start_fwd;#0 based for array
				my @arr_lig_probe_sequence =  split(//, $lig_probe_sequence_fwd);
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$lig_probe_sequence_fwd = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_sequence_fwd = $lig_probe_sequence_fwd . $bases;
				}
				my $allele_probe_position = $allele_probe_seq_position + 1;
				print OUTFILE "$mip_count\t$chr\t$ext_probe_start_fwd\t$ext_probe_stop_fwd\t$ext_probe_sequence_fwd\t$lig_probe_start_fwd\t$lig_probe_stop_fwd\t$lig_probe_sequence_fwd\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
				print "$line_count/$number_of_features\t$mip_count\n";
			}
		}
	}
	else
	{
		$bad_resriction_site_design_count++;
	}

	#check reverse for restriction sites
	$ntalwl_seq_5_count = 0;
	$ntalwl_seq_3_count = 0;
	$nbbsrdi_seq_5_count = 0;
	$nbbsrdi_seq_3_count = 0;
	$ext_probe_sequence_rev = uc(substr($chr_sequence,$ext_probe_start_rev - 1,$extension_arm_length));
	$lig_probe_sequence_rev = uc(substr($chr_sequence,$lig_probe_start_rev - 1,$ligation_arm_length));
	&reverse_complement_target;
	&reverse_complement_arms;
	$scan_target_sequence_with_probe_arms = $ext_probe_sequence_rev . $scan_target_sequence . $lig_probe_sequence_rev;
	#print "$scan_target_sequence_with_probe_arms\n";
	#CCTCCCCTGGGCCCCCCACCGGCACCCTCCGCCGCCCCTTCTTGAACACACTCAGGAATGGCTGCTCCAGGTTTTCTCGCTGGTTCTCCAGGTCCAGGATGTCCTAGGAGGAGTAGAGCTCAGGGGAGGGGGCTTTTCCCAGGTTCCTCACA
	$mip_seq = $universal_left_end_mip_seq . $lig_probe_sequence_rev . $universal_middle_mip_seq . $ext_probe_sequence_rev . $universal_right_end_mip_seq;
	#print "$mip_seq\n";
	#AGGACCGGATCAACTCTTTTCCCAGGTTCCTCACACTTCAGCTTCCCGATATCCGACGGTAGTGTCCTCCCCTGGGCCCCCCACCCATTGCGTGAACCGA
	#check the mip for restriction sites
	$length_mip_seq = length($mip_seq) - 1;#0 based
	for my $int (0 .. ($length_mip_seq - 4))
	{
		my $sub_mip_seq_length_5 = substr($mip_seq, $int, 5);
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_mip_seq - 5))
	{
		my $sub_mip_seq_length_6 = substr($mip_seq, $int, 6);
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}

	#check the target with probe arms
	my $length_target_sequence_with_probe_arms = length($scan_target_sequence_with_probe_arms) - 1;#0 based
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 4))
	{
		my $sub_target_seq_length_5 = substr($length_target_sequence_with_probe_arms, $int, 5);
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 5))
	{
		my $sub_target_seq_length_6 = substr($length_target_sequence_with_probe_arms, $int, 6);
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}

	if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
	{
		#do the rev/- strand snps
		my $snp_count_rev = 0;
		for my $int ($ext_probe_start_rev .. $ext_probe_stop_rev)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count_rev++;
				$snp_position = $int;
			}
		}
		for my $int ($lig_probe_start_rev .. $lig_probe_stop_rev)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count_rev++;
				$snp_position = $int;
			}
		}
		if ($snp_count_rev == 0)
		{
			$mip_count++;
			print OUTFILE "$mip_count\t$chr\t$ext_probe_start_rev\t$ext_probe_stop_rev\t$ext_probe_sequence_rev\t$lig_probe_start_rev\t$lig_probe_stop_rev\t$lig_probe_sequence_rev\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\n";
			print "$line_count/$number_of_features\t$mip_count\n";
		}
		elsif ($snp_count_rev == 1)
		{#make for allele_2
			#make two mips, one on each strand for each allele of the snp single snp found
			&get_alleles_for_position;#$snp_position
			if ($snp_position >= $ext_probe_start_rev && $snp_position <= $ext_probe_stop_rev)
			{
				&reverse_complement_arms;#put it back in genome orientation for the snp change
				$mip_count++;
				my $allele_probe_seq_position = $snp_position - $ext_probe_start_rev;#0 based for array
				my @arr_ext_probe_sequence =  split(//, $ext_probe_sequence_rev);
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$ext_probe_sequence_rev = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_sequence_rev = $ext_probe_sequence_rev . $bases;
				}
				&reverse_complement_arms;
				my $allele_probe_position = $extension_arm_length - ($allele_probe_seq_position + 1) + 1;
				print OUTFILE "$mip_count\t$chr\t$ext_probe_start_rev\t$ext_probe_stop_rev\t$ext_probe_sequence_rev\t$lig_probe_start_rev\t$lig_probe_stop_rev\t$lig_probe_sequence_rev\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
				print "$line_count/$number_of_features\t$mip_count\n";
			}
			elsif ($snp_position >= $lig_probe_start_rev && $snp_position <= $lig_probe_stop_rev)						
			{#must be between the ligerse target arms
				&reverse_complement_arms;#put it back in genome orientation for the snp change
				$mip_count++;
				my $allele_probe_seq_position = $snp_position - $lig_probe_start_rev;#0 based for array
				my @arr_lig_probe_sequence =  split(//, $lig_probe_sequence_rev);
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$lig_probe_sequence_rev = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_sequence_rev = $lig_probe_sequence_rev . $bases;
				}
				&reverse_complement_arms;
				my $allele_probe_position = $ligation_arm_length - ($allele_probe_seq_position + 1) + 1;
				print OUTFILE "$mip_count\t$chr\t$ext_probe_start_rev\t$ext_probe_stop_rev\t$ext_probe_sequence_rev\t$lig_probe_start_rev\t$lig_probe_stop_rev\t$lig_probe_sequence_rev\t$mip_scan_start_position\t$mip_scan_stop_position\t$scan_target_sequence\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
				print "$line_count/$number_of_features\t$mip_count\n";
			}
		}
	}
	else
	{
		$bad_resriction_site_design_count++;
	}
}
sub reverse_complement_arms
{
	#do the ext
	my $flipped_seq;
	my $x;
	while ($ext_probe_sequence_rev ne "")
	{
		$x = chop($ext_probe_sequence_rev);
		$flipped_seq .= $x;
	}
	### complement sequence ###
	$flipped_seq =~ tr/ATCG/TAGC/;
	$ext_probe_sequence_rev = $flipped_seq;
	#do the lig
	$flipped_seq = "";
	$x = "";
	while ($lig_probe_sequence_rev ne "")
	{
		$x = chop($lig_probe_sequence_rev);
		$flipped_seq .= $x;
	}
	### complement sequence ###
	$flipped_seq =~ tr/ATCG/TAGC/;
	$lig_probe_sequence_rev = $flipped_seq;
}

sub reverse_complement_target
{
	my $flipped_seq;
	my $x;
	while ($scan_target_sequence ne "")
	{
		$x = chop($scan_target_sequence);
		$flipped_seq .= $x;
	}
	$flipped_seq =~ tr/ATCG/TAGC/;
	$scan_target_sequence = $flipped_seq;
}

sub get_alleles_for_position
{
	my @alleles =  split(//, $chr_snp_positions{$chr}{$snp_position});
	$allele_1 = $alleles[0];
	$allele_2 = $alleles[1];
}
sub load_chr_snps
{
	my $snp_load_count = 0;
	#my $snp_file = "/nfs/home/bard/dbSNP/build_129/xml_chr/gt_chr".$chr.".xml_trimmed_to_site_info.rs_orientation_info";
	#my $snp_file = "/nfs/home/bard/dbSNP/build_129/ds_xml_chr/ds_ch".$chr.".xml";
	
	############## COMMENTED OUT BY JOE ON 110613
	# my $snp_file = "/net/grc/vol1/references/human/snps/dbSNP/GATK_rod_hg18_build_130/dbsnp_130_hg18.rod";
	
	open (SNPFILE, $snp_file) or die "Can't open\n $snp_file\n";
	while (<SNPFILE>)
	{
		my $line = $_;
		chomp $line;
		my @line_contents = split(/\s+/, $line);
		my $chr = $line_contents[1];
		$chr =~ s/chr//;
		my $position = $line_contents[2]+1; #BJO
		my $snp_strand_orientation = $line_contents[6];
		my $alleles = $line_contents[9];
		my @alleles_content = split(/\//, $alleles);
		my $a1 = $alleles_content[0];
		my $a2 = $alleles_content[1];
		if ($snp_strand_orientation eq "-")
		{#flip them to be forward
			$a1 =~ tr/ATCG/TAGC/;
			$a2 =~ tr/ATCG/TAGC/;
			$chr_snp_positions{$chr}{$position} = $a1 . $a2;
			#print "$chr\t$position\t$a1\t$a2\tlig->ext\n";
			$snp_load_count++;
			$total_snp_load_count++;
			#print "$total_snp_load_count\n";
		}
		else#just load the forward strand
		{
			$chr_snp_positions{$chr}{$position} = $a1 . $a2;
			#print "$chr\t$position\t$a1\t$a2\text\n";
			$snp_load_count++;
			$total_snp_load_count++;
			#print "$total_snp_load_count\n";
		}
	}
	close SNPFILE;
}
sub parseCommandLine 
{
    my ( $useage ) = "

-regions_to_scan

format:
#bed file
#chromosome start(0-based) stop(1-based)
chr1 58952 59873
chr1 357520 358462
chr1 610957 611899
chr1 851183 851258
chr1 855396 855581
chr1 856280 856334

-mip_scan_size 112

(2 x 76 bp reads on a GAII) - (2 X target_arm_length)

-extension_arm_length 19

	generally is less than the ligation_arm_length

-ligation_arm_length 21

-feature_flank 56

	optional, defaults to 0

-check_for_snps y/n

-genome_dir /net/shendure/...

-snp_file /net/shendure/...

";

        for (my $i = 0; $i <= $#ARGV; $i++)
        {
                if ($ARGV[$i] =~ /^-/)
                {
                        $arg{$ARGV[$i]} = $ARGV[$i+1];
                }
        }
        die($useage) if (!($arg{-regions_to_scan}));
        die($useage) if (!($arg{-ligation_arm_length}));
        die($useage) if (!($arg{-extension_arm_length}));
        die($useage) if (!($arg{-mip_scan_size}));
        die($useage) if (!($arg{-check_for_snps}));
		die($useage) if (!($arg{-genome_dir}));
		die($useage) if (!($arg{-snp_file}));
}
sub get_regions_to_scan
{
	open (FILE, $regions_to_scan);
	@regions_to_scan_contents = <FILE>;
	close FILE;
}
sub get_chr_fasta_sequence
{
	open (FILE, $chr_fasta) or die "Can't open $chr_fasta\n";
	my @contents = <FILE>;
	close FILE;
	$chr_sequence = "";
	foreach my $lines (@contents)
	{
		chomp $lines;
		if ($lines =~ m/>/g)
		{}
		else
		{
			$chr_sequence = $chr_sequence . uc($lines);
		}
	}
}
