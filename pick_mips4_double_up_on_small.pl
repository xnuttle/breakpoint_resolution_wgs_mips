#!/usr/bin/perl
######################################################
#written by Ian Byrell Stanaway bard@u.washington.edu#
######################################################
use strict;
use warnings;
use POSIX qw(ceil floor);
my %arg;
&parseCommandLine;
my $bad_junctions_file = $arg{-bad_junctions_file};
my %bad_junctions;
&get_bad_junctions_file;
my $regions_to_scan = $arg{-regions_to_scan};
my @regions_to_scan_contents;
&get_regions_to_scan;
my $mips_file_list = $arg{-mips_file_list};
my %chr_target_start_mip;
my $mips_quality_file_list = $arg{-mips_quality_file_list};
my $file_id_to_get;
my %mip_index_rank;
my $mip_overlap = $arg{-mip_overlap};
my $mip_scan_size = $arg{-mip_scan_size};
my $feature_flank = 0;
if (exists $arg{-feature_flank})
{
	$feature_flank = $arg{-feature_flank};
}

&get_mips_files;

my $feature_start_position;
my $feature_stop_position;

my $seq_arm_5;
my $seq_arm_3;
my %chr_pos_strand_used_arm_bases;

my $feature_mip_count = 0;
my $probe_copy_number = 0;
my $picked_mip_count = 0;
my $mip_strand;
my $mip_type;# t tiled or s single
my $mip_search_space;#bp we can move/wiggle and still have the target in the mip target zone

my $mip_target_start_position;
my $mip_target_stop_position;
my $last_mip_target_stop_position = "NA";

my $snp_position;
my $allele_1;
my $allele_2;

my $mip_finished;# y/n
my $design_attempts = 0;
my $bad_design_count = 0;
my $design_done = 0;
my $chr;

#my $outfile = "$regions_to_scan" . ".picked_mip_probe_arms.$chromosome_to_scan.$mips_file_list";
$regions_to_scan =~ s/\r//g;
$mips_file_list =~ s/\r//g;

my $mip_basename = $mips_file_list;

if($mips_file_list =~ m/^\S+\/\S+/)
{
	$mips_file_list =~ m/.+\//;
	$mip_basename = $';
}

my $outfile = "$regions_to_scan" . ".picked_mip_probe_arms.$mip_basename";
open (OUTFILE, ">$outfile") or die "Cannot Open $outfile";
print OUTFILE ">mip_pick_count\trank_score\tchr\text_probe_start\text_probe_stop\text_probe_sequence\text_copy_count\tlig_probe_start\tlig_probe_stop\tlig_probe_sequence\tlig_copy_count\tmip_target_start_position\tmip_target_stop_position\tmip_target_sequence\tfeature_start_position\tfeature_stop_position\tfeature_mip_count\tprobe_strand\tnotes\n";
my $line_count = 0;
my $number_of_features = $#regions_to_scan_contents;
#print "picking mips\n";
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
		$chr =~ s/chr//gi;
		$last_mip_target_stop_position = "NA";
		$feature_start_position = $line_contents[1] - $feature_flank + 1;#make 1 based, is 0 based in file, #BJO
		$feature_stop_position = $line_contents[2] + $feature_flank;
		my $feature_length = $feature_stop_position - $feature_start_position + 1;
		my $number_of_mips = int($feature_length / $mip_scan_size) + 1;
		#print "large $feature_length bp feature, designing ~ $number_of_mips mips\n";
		$mip_target_start_position = 0;
		$mip_target_stop_position = 0;
		$feature_mip_count = 0;
		$mip_strand = "+";# can be +/- and alternates every other mip
		if ($last_mip_target_stop_position eq "NA")
		{
			$last_mip_target_stop_position = $feature_start_position;
		}
		#while ($mip_target_stop_position <= $feature_stop_position)
		while ($last_mip_target_stop_position < $feature_stop_position)
		{
			$feature_mip_count++;
			if ($feature_mip_count > 1)
			{#get the last probe position
				$mip_target_start_position = $last_mip_target_stop_position - $mip_overlap + 1;
				if ($mip_strand eq "+")#strand alternation
				{
					$mip_strand = "-";
				}
				else
				{
					$mip_strand = "+";
				}
			}
			else
			{
				$mip_target_start_position = $feature_start_position;
			}
			$mip_target_stop_position = $mip_target_start_position + $mip_scan_size - 1;
			$mip_search_space = $mip_scan_size / 2;
			&pick_mip_in_search_space;
		}

		if ($feature_mip_count == 1)#make another so we have 2 mips, one on each strand
		{
			$last_mip_target_stop_position = $feature_start_position;
			$mip_target_start_position = $feature_start_position;
			if ($mip_strand eq "+")#put it on the other strand
			{
				$mip_strand = "-";
			}
			else
			{
				$mip_strand = "+";
			}
			$mip_target_stop_position = $mip_target_start_position + $mip_scan_size - 1;
			$mip_search_space = $mip_scan_size / 2;
			&pick_mip_in_search_space;
		}
	}
}
close OUTFILE;

sub pick_mip_in_search_space
{
	my $snp_count = 0;
	my $design_done = 0;
	my $best_mip = "NA";
	my $best_target_start_position;
	my $best_target_stop_position;
	my $best_rank = -100000;#-1 is lowest in file
	
	my $best_snp_mip_strand_1 = "NA";
	#my $best_snp_mip_strand_2 = "NA";
	my $best_snp_target_start_position;
	my $best_snp_target_stop_position;
	my $best_snp_rank = -100000;#-1 is lowest in file

	my $tmp_mip_target_start_position = $mip_target_start_position;
	my $first_mip = "NA";
	my $first_mip_rank;
	my $counter = 0;
	while ($design_done == 0)
	{
		$counter++;
		my $mip = "NA";
		my $rank = -200000000;
		my @mip_contents;
		if (exists $chr_target_start_mip{$chr}{$tmp_mip_target_start_position}{$mip_strand})
		{
			$mip = $chr_target_start_mip{$chr}{$tmp_mip_target_start_position}{$mip_strand};
			@mip_contents = split(/\s+/, $mip);
			$rank = $mip_contents[1];
		}
		my $used_arm_space_flag = 0;
		if ($mip ne "NA")
		{
			my $chr = $mip_contents[2];		
			my $ext_probe_start = $mip_contents[3];
			my $ext_probe_stop = $mip_contents[4];
			my $lig_probe_start = $mip_contents[7];
			my $lig_probe_stop = $mip_contents[8];
			my $strand = $mip_contents[17];

			foreach my $pos ($ext_probe_start .. $ext_probe_stop)
			{
				#print "$mip\n$chr\t$pos\t$mip_strand\n";
				if (exists $chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand})
				{
					$used_arm_space_flag = 1;
				}
			}
			foreach my $pos ($lig_probe_start .. $lig_probe_stop)
			{
				#print "$chr\t$pos\t$mip_strand\n";
				if (exists $chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand})
				{
					$used_arm_space_flag = 1;				
				}
			}
		}
		if ($first_mip eq "NA" && $mip !~ m/snp/ && $mip ne "NA" && $used_arm_space_flag != 1)
		{
			$first_mip = $mip;
			$first_mip_rank = $rank;
		}

		if ($mip eq "NA" || $rank == -200000000 || $used_arm_space_flag == 1)
		{
			#if ($used_arm_space_flag == 1)
			#{
			#	print "invalidating arm\n";
			#}
		}
		else
		{
			my $ext_arm_copy_count = $mip_contents[6];
			my $lig_arm_copy_count = $mip_contents[10];
			my $copy_count;
			if ($lig_arm_copy_count > $ext_arm_copy_count)
			{
				$copy_count = $lig_arm_copy_count;
			}
			elsif ($lig_arm_copy_count < $ext_arm_copy_count)
			{
				$copy_count = $ext_arm_copy_count;
			}
			else
			{
				$copy_count = $ext_arm_copy_count;
			}
			if ($copy_count >= 10)
			{
				$rank = $rank - 10000;
			}
			if ($rank > $best_rank)
			{
				if ($mip =~ m/snp/)
				{}
				else
				{
					$best_mip = $mip;
					$best_rank = $rank;
					$best_target_start_position = $tmp_mip_target_start_position;
					$best_target_stop_position = $mip_contents[12];
				}
			}
			if ($rank > $best_snp_rank)
			{
				if ($mip =~ m/snp/)
				{
					$best_snp_mip_strand_1 = $mip;
					$best_snp_rank = $rank;
					$best_snp_target_start_position = $tmp_mip_target_start_position;
					$best_snp_target_stop_position = $mip_contents[12];
				}
			}
		}
		$tmp_mip_target_start_position--;
		$mip_search_space--;
		if ($mip_search_space <= 0)
		{
			$design_done = 1;
		}
	}
	if ($best_mip ne "NA")
	{
		$picked_mip_count++;
		$last_mip_target_stop_position = $best_target_stop_position;
		$best_mip =~ s/^\d+/$picked_mip_count/g;
		print OUTFILE "$best_mip\tbest_mip\n";
		#print "$line_count/$number_of_features\t$picked_mip_count\n";
		my @mip_contents = split(/\s+/, $best_mip);
		my $chr = $mip_contents[2];		
		my $ext_probe_start = $mip_contents[3];
		my $ext_probe_stop = $mip_contents[4];
		my $lig_probe_start = $mip_contents[7];
		my $lig_probe_stop = $mip_contents[8];
		my $strand = $mip_contents[17];
		foreach my $pos ($ext_probe_start .. $ext_probe_stop)
		{
			$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
		}
		foreach my $pos ($lig_probe_start .. $lig_probe_stop)
		{
			$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
		}
	}
	#elsif ($best_snp_mip_strand_1 ne "NA" && $best_snp_mip_strand_2 ne "NA" && $best_snp_mip_strand_1 ne "" && $best_snp_mip_strand_2 ne "")
	elsif ($best_snp_mip_strand_1 ne "NA")
	{
		$picked_mip_count++;
		$last_mip_target_stop_position = $best_snp_target_stop_position;
		$best_snp_mip_strand_1 =~ s/^\d+/$picked_mip_count/g;
		print OUTFILE "$best_snp_mip_strand_1:picked_snp\n";
		my @mip_contents = split(/\s+/, $best_snp_mip_strand_1);
		my $chr = $mip_contents[2];		
		my $ext_probe_start = $mip_contents[3];
		my $ext_probe_stop = $mip_contents[4];
		my $lig_probe_start = $mip_contents[7];
		my $lig_probe_stop = $mip_contents[8];
		my $strand = $mip_contents[17];
		foreach my $pos ($ext_probe_start .. $ext_probe_stop)
		{
			$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
		}
		foreach my $pos ($lig_probe_start .. $lig_probe_stop)
		{
			$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
		}
	}
	else
	{
		if ($first_mip eq "NA")
		{
			my $mip_target_start_position_missing = $last_mip_target_stop_position + 1;
			$last_mip_target_stop_position = $last_mip_target_stop_position + ($mip_scan_size / 2);#move half a mip since we didn't test the starts, and we can find anything in the first half (search space)
			$picked_mip_count++;
			#print OUTFILE "$picked_mip_count\tNA\t$chromosome_to_scan\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$mip_target_start_position_missing\t$last_mip_target_stop_position\tNA\tNA\tno_mips_found_in_search_space\n";
			print OUTFILE "$picked_mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$mip_target_start_position_missing\t$last_mip_target_stop_position\tNA\tNA\tno_mips_found_in_search_space\n";
			#print "$line_count/$number_of_features\t$picked_mip_count\n";
		}
		else
		{
			my @first_mip_contents =  split(/\s+/, $first_mip);
			my $first_target_stop_position = $first_mip_contents[12];
			$last_mip_target_stop_position = $first_target_stop_position;
			$picked_mip_count++;
			$first_mip =~ s/^\d+/$picked_mip_count/g;
			print OUTFILE "$first_mip\tfirst_mip_used\n";
			#print "$line_count/$number_of_features\t$picked_mip_count\n";
			my @mip_contents = split(/\s+/, $first_mip);
			my $chr = $mip_contents[2];		
			my $ext_probe_start = $mip_contents[3];
			my $ext_probe_stop = $mip_contents[4];
			my $lig_probe_start = $mip_contents[7];
			my $lig_probe_stop = $mip_contents[8];
			my $strand = $mip_contents[17];
			foreach my $pos ($ext_probe_start .. $ext_probe_stop)
			{
				$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
			}
			foreach my $pos ($lig_probe_start .. $lig_probe_stop)
			{
				$chr_pos_strand_used_arm_bases{$chr}{$pos}{$strand} = 0;
			}
		}
	}
}

sub get_mips_files
{
	open (FILE, $mips_file_list) or die "can't find the mip file list";
	my @mips_file_list_contents = <FILE>;
	close FILE;
	my $number_of_mip_files = $#mips_file_list_contents + 1;
	my $mip_file_count = 0;
	foreach my $mips_file (@mips_file_list_contents)
	{
		$mip_file_count++;
		chomp $mips_file;
		my $mip_line_count = 0;
		print "$mip_file_count / $number_of_mip_files loading mips $mips_file \n";
		open (FILE, $mips_file);
		while (<FILE>)
		{
			my $lines = $_;
			if ($lines =~ m/^>|^#|^\n|^\s/g)#header
			{}
			else
			{
				chomp $lines;
				$mip_line_count++;
				#print $mip_line_count . " $mips_file\n";
				my @line_contents = split(/\s+/, $lines);
				my $chromosome = $line_contents[2];
				my $target_start_position = $line_contents[11];
				my $strand = $line_contents[17];
				my $arm_ext = $line_contents[5];
				my $arm_lig = $line_contents[9];
				my $rank = $line_contents[1];
				my $gap_fill = $line_contents[13];
				my $arm_lig_2mer =  substr($arm_lig, 0, 2);
				my $gap_fill_2mer =  substr($gap_fill, 110, 2);
				my $junction = $gap_fill_2mer . $arm_lig_2mer;
				if (exists $bad_junctions{$junction})
				{#T junctions are bad for the ligation reaction
				}
				else
				{
					if ($arm_ext =~ m/-|n/gi || $arm_lig =~ m/-|n/gi || $rank eq "" || $lines =~ m/snp:\d+:\d+:-\/\w+|snp:\d+:\d+:\w+\/-/)
					{}
					elsif ($lines =~ m/snp/)
					{
						if (exists $chr_target_start_mip{$chromosome}{$target_start_position}{$strand})
						{
							my $old_mip = $chr_target_start_mip{$chromosome}{$target_start_position}{$strand};
							my @old_mip_contents = split(/\s+/, $old_mip);
							my $old_rank = $old_mip_contents[1];
							if ($rank > $old_rank && $old_mip =~ m/snp/)
							{
								#print "replacing mip";
								$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $lines;
							}
							elsif ($rank == $old_rank && $old_mip =~ m/snp/)
							{
								$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $chr_target_start_mip{$chromosome}{$target_start_position}{$strand} . "|" . $lines;
							}
						}
						else
						{
							$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $lines;
						}
					}
					else
					{
						if (exists $chr_target_start_mip{$chromosome}{$target_start_position}{$strand})
						{
							my $old_mip = $chr_target_start_mip{$chromosome}{$target_start_position}{$strand};
							my @old_mip_contents = split(/\s+/, $old_mip);
							my $old_rank = $old_mip_contents[1];
							if ($rank > $old_rank)
							{
								$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $lines;
							}
							elsif ($rank == $old_rank)
							{
								$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $chr_target_start_mip{$chromosome}{$target_start_position}{$strand} . "|" . $lines;
							}
							if ($old_mip =~ m/snp/)
							{
								#print "replacing mip";
								$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $lines;
							}
						}
						else
						{
							$chr_target_start_mip{$chromosome}{$target_start_position}{$strand} = $lines;
						}
					}
				}
			}
		}
		close FILE;
		#print "$mips_file_list loaded mips\n";
	}

	foreach my $chromosomes (keys %chr_target_start_mip)
	{#and randomly pick one of them
		foreach my $target_start_position (keys %{$chr_target_start_mip{$chromosomes}})
		{
			foreach my $strand (keys %{$chr_target_start_mip{$chromosomes}{$target_start_position}})
			{
				if ($chr_target_start_mip{$chromosomes}{$target_start_position}{$strand} =~ m/\|/)#is more than one mip
				{
					#print "has pipe $chr_target_start_mip{$chromosomes}{$target_start_position}{$strand}\n";
					my $mips = $chr_target_start_mip{$chromosomes}{$target_start_position}{$strand};
					#print "has pipe2 $mips\n";
					my @possible_mips = split(/\|/, $mips);
					my $possible_mip_count = $#possible_mips + 1;
					my $random_number = rand($possible_mip_count);#between 0 - $possible_mip_count
					my $ceil_random_number = ceil($random_number) - 1;#round it up to the integer
					my $mip_to_use = $possible_mips[$ceil_random_number];
					#print "array pos $ceil_random_number  $random_number $possible_mip_count  mip_to_use\t$mip_to_use\n";
					$chr_target_start_mip{$chromosomes}{$target_start_position}{$strand} = $mip_to_use;
				}
			}
		}
	}
}

sub parseCommandLine 
{
    my ( $useage ) = "

#-chromosome_to_scan

-regions_to_scan

format:
#coding regions plus 2 bases each side for CCDS.20080902.txt, unique regions concatenated
#chromosome start(0-based) stop(1-based)
1 58952 59873
1 357520 358462
1 610957 611899
1 851183 851258
1 855396 855581
1 856280 856334

-mips_file_list

-mip_scan_size 112

(2 x 76 bp reads on a GAII) - (2 X target_arm_length)

-mip_overlap 1

-feature_flank 2

	optional, defaults to 0

-bad_junctions_file

";

        for (my $i = 0; $i <= $#ARGV; $i++)
        {
                if ($ARGV[$i] =~ /^-/)
                {
                        $arg{$ARGV[$i]} = $ARGV[$i+1];
                }
        }
        die($useage) if (!($arg{-bad_junctions_file}));
        die($useage) if (!($arg{-regions_to_scan}));
        die($useage) if (!($arg{-mips_file_list}));
        die($useage) if (!($arg{-mip_overlap}));
        die($useage) if (!($arg{-mip_scan_size}));
}
sub get_regions_to_scan
{
	open (FILE, $regions_to_scan);
	@regions_to_scan_contents = <FILE>;
	close FILE;
}
sub get_bad_junctions_file
{
	open (FILE, $bad_junctions_file) or die "can't find the mip file list";
	while (<FILE>)
	{
		my $lines = $_;
		if ($lines =~ m/^>|^#|^\n|^\s/g)#header
		{}
		else
		{
			chomp $lines;
			#print $mip_line_count . " $lines\n";
			my @line_contents = split(/\s+/, $lines);
			my $junction = $line_contents[4];
			$bad_junctions{$junction} = 0;
		}
	}
	close FILE;
}
