//Xander Nuttle
//tagged_mrfast_output_to_mipcounts.c
//Prints to standard output paralog-specific read counts for each MIP target for a single individual in an experiment
//Call: ./tagged_mrfast_output_to_mipcounts miptargets_file text_file_with_names_of_mapping_output_files_for_an_individual (int)molecular_tag_length (int)individual_number
//
//This program takes in (1) a file containing information on MIP targets, (2) a file with names of all gzipped mapping output 
//SAM files from mrFAST for a particular individual, (3) the molecular tag length, and (4) the number of the individual in the 
//experiment (first individual to be analyzed gets the number 1) and outputs to standard output paralog-specific read counts
//for each MIP target for each individual, where the final counts are informed by observing distinct molecular tags.
//
//For each read pair, the program performs the following steps to assess mapping and increment the appropriate read count:
//	-determine if mapping location corresponds to a valid MIP target (if not, go to next read pair)
//	-determine if insert size calculated from mapping is concordant with expectation from MIP design (if not, go to next read pair)
//	-ensure base quality is >= Phred 30 at all bases mapped to sites having variation between paralogs in contig target sequences (if not, go to next read pair)
//	-ensure no mismatches at bases mapped to likely fixed (high-scoring) SUN sites (if not, go to next read pair)
//	-determine paralog-specificity for MIP target the read pair mapped to
//	-determine whether the molecular tag has already been seen for this MIP-paralog target, and if it hasn't, increment the appropriate MIP_target-paralog read count
//
//After all mapped read pairs have been processed, the program prints the counts to standard output.
//
//THIS PROGRAM DOES NOT ANALYZE SEQUENCE CONTENT!!!
//
//This program should be run from the directory containing mrfast output files (gzipped sam files), i.e. a "mrfast_mapping_output" directory

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<zlib.h>
#define KYEL "\x1B[33m"
#define KRED "\x1B[31m"

int parse_cigar_and_md(char*cigar_string,char*md_string,char*read_sequence,char*original_quality_array,long mapping_location,long expected_mapping_location,char*sequence_array,char*quality_array,int*is_base_mm_array);

int main(int argc,char*argv[])
{
	//check to make sure there are there are the correct number of command line arguments provided
	if(argc!=5) //the program call counts as the 1st command line argument; thus, with 4 necessary inputs, argc must be 5
	{
		printf(KRED "ERROR due to improper input!!! The program call should take the following format:\n");
		printf(KYEL "PATH_TO_PROGRAM/tagged_mrfast_output_to_mipcounts miptargets_file text_file_with_names_of_mapping_output_files_for_an_individual (int)molecular_tag_length (int)individual_number\n");
		printf(KYEL "\nProgram Help:\n\n");
		system("head -23 /net/gs/vol2/home/xnuttle/MIPs/analysis_programs/tagged_mrfast_output_to_mipcounts.c");
		return 1;
	}

	//get information about number of mip targets and maximum number of contigs per family
	long max_num_contigs=0,num_mip_targets=0;
	FILE*miptargetsfile;
	fpos_t pos;
	char dummy[501];
	char specstring[51];	
	miptargetsfile=fopen(*(argv+1),"r");
  fgetpos(miptargetsfile,&pos);
	while(fscanf(miptargetsfile,"%s %s %s %s %s %s %s %s %s",dummy,dummy,dummy,dummy,dummy,dummy,specstring,dummy,dummy)==9)
  {
		num_mip_targets++;
		if(strlen(specstring)>max_num_contigs)
		{
			max_num_contigs=strlen(specstring);
		}
	}	
	fsetpos(miptargetsfile,&pos);
	
	//set up mip structure
	struct mip
	{
		char master_target_name[51];
		long master_coordinate;
		char contig_names[max_num_contigs][51];
		long contig_start_coords[max_num_contigs];
		char mip_type;
		int important_bases[152];
		int specificity[max_num_contigs];
		char orientation;
	};

	//set up structure to store information on molecular tags seen and mip paralog-specific sequence read counts
	long taglength=strtol(*(argv+3),NULL,10);
	struct mipdata
	{
		long counts[num_mip_targets][max_num_contigs+1];
		char tags[num_mip_targets][max_num_contigs+1][5000][taglength+1]; //highly unlikely more than 5000 distinct tags will be observed at a single MIP target
	};

	//allocate memory to store mip information
	struct mip*mips;
	mips=(struct mip*)malloc(num_mip_targets*sizeof(struct mip));
	if(mips==NULL)
	{
		printf("Memory allocation for mip information failed!\n");
		return 1;
	}

	//allocate memory to store data for paralog-specific read counts and molecular tags seen
	struct mipdata*data;
	data=(struct mipdata*)malloc(sizeof(struct mipdata));
	if(data==NULL)
	{
		printf("Memory allocation for mip read count and molecular tag data failed!\n");
    return 1;
	}

	//read in file containing information on mip target locations and store data in mip structure elements
	char line[5001];
	char ch;
	int j,y,z;
	long master_start,master_end,var_position;
	for(j=0;j<num_mip_targets;j++)
	{
		//get first field - master sequence name - and store in mip structure
		z=0;
		while((ch=getc(miptargetsfile))!='\t')
		{
			mips[j].master_target_name[z]=ch;
			z++;
		}
		mips[j].master_target_name[z]='\0';
		//get second field - master sequence start and end coordinates - and store an average of these in mip structure
		fscanf(miptargetsfile,"%ld%c%ld%c",&master_start,&ch,&master_end,&ch);
		mips[j].master_coordinate=(master_start+master_end)/2;
		//get third field - contig names and start coordinates - and store in mip structure
		for(y=0;y<max_num_contigs;y++)
		{
			strcpy(mips[j].contig_names[y],"");
			mips[j].contig_start_coords[y]=0;
		}
		y=0;
		while((ch=getc(miptargetsfile))!='\t')
		{
			z=0;
			while(ch!=':')
			{
				mips[j].contig_names[y][z]=ch;
				z++;
				ch=getc(miptargetsfile);
			}
			mips[j].contig_names[y][z]='\0';
			z=0;
			ch=getc(miptargetsfile);
			while(ch!=';')
			{
				line[z]=ch;
				z++;
				ch=getc(miptargetsfile);
			}
			line[z]='\0';
			mips[j].contig_start_coords[y]=strtol(line,NULL,10);
			y++;
		}
		//get fourth field - mip type - and store in mip structure
		mips[j].mip_type=getc(miptargetsfile);
		ch=getc(miptargetsfile);
		//get 5th and 6th fields - local variant base and fixed SUN locations - and store in mip structure in a 152-mer sequence of 0s,1s, and numbers >1
		//0=base invariant at this position,1=some variation between paralogs at this position but not a fixed SUN,2+=fixed SUN at this position
		for(z=0;z<152;z++)
		{
			mips[j].important_bases[z]=0;
		}
		while((ch=getc(miptargetsfile))!='\t')
		{
			z=0;
			while(ch!=',')
			{
				line[z]=ch;
				z++;
				ch=getc(miptargetsfile);
			}
			line[z]='\0';
			var_position=strtol(line,NULL,10);
			if(var_position!=0)
			{
				mips[j].important_bases[(var_position-1)]++;
			}
		}
		while((ch=getc(miptargetsfile))!='\t')
    {
      z=0;
      while(ch!=',')
      {
        line[z]=ch;
        z++;
        ch=getc(miptargetsfile);
      }
      line[z]='\0';
      var_position=strtol(line,NULL,10);
      if(var_position!=0)
			{
				mips[j].important_bases[var_position-1]++;
			}
    }
		//get 7th and 8th fields - specificity and orientation - and store in mip structure
		for(z=0;z<max_num_contigs;z++)
		{
			mips[j].specificity[z]=0;
		}
		z=0;
		while((ch=getc(miptargetsfile))!='\t')
		{
			line[0]=ch;
			line[1]='\0';
			var_position=strtol(line,NULL,10);
			mips[j].specificity[z]=(int)var_position;
			z++;
		}
		mips[j].orientation=getc(miptargetsfile);
		//ignore 9th field - number of bases in 1st MIP arm
		while(ch=getc(miptargetsfile)!='\n')
			continue;
	}
	fclose(miptargetsfile);
	
	//initialize mip paralog-specific read counts and molecular tag information
	int k;
	for(j=0;j<num_mip_targets;j++)
	{
		for(k=0;k<(max_num_contigs+1);k++)
		{
			data->counts[j][k]=0;
			for(y=0;y<5000;y++)
      {
        for(z=0;z<taglength;z++)
        {
          data->tags[j][k][y][z]='N';
        }
        data->tags[j][k][y][z]='\0';
      }
		}
	}

	//read in text file containing names of all gzipped mapping output files to input to program, and process each file one-by-one
	FILE*textfile;
	int i;
	char input_file_name[51];
	char line2[5001];
  char mapped_orientation;
	char tagseq[taglength+1];
	tagseq[taglength]='\0';
  char mapped_contig[18];
  long mapping_loc,mapping_loc2;
  long target_size;
  char quality[152],original_quality[76],original_quality2[76];
	char sequence[152],read1seq[76],read2seq[76];
  int is_base_mm[152];
  char cigar[301],cigar2[301];
  char md[301],md2[301];
  char*tab_locations1[12];
  char*tab_locations2[12];
	textfile=fopen(*(argv+2),"r");
	while(fscanf(textfile,"%s",input_file_name)==1)
	{
		//process gzipped file containing mapped reads
		gzFile*in;
		in=gzopen(input_file_name,"r");
		while((ch=gzgetc(in))!=EOF)
		{
			i=0;
			line[0]=ch;
			z=1;
			while((ch=gzgetc(in))!='\n')
			{
				line[z]=ch;
				if(ch=='\t')
				{
					tab_locations1[i]=line+z;
					i++;
				}
				z++;
			}
			line[z]='\0';
			y=0;
			i=0;
			while((ch=gzgetc(in))!='\n')
  		{
    		line2[y]=ch;
    		if(ch=='\t')
    		{
      		tab_locations2[i]=line2+y;
      		i++;
    		}
    		y++;
  		}
			line2[y]='\0';
			strncpy(tagseq,strchr(line,'$')+1,taglength);
			strncpy(mapped_contig,tab_locations1[1]+1,tab_locations1[2]-tab_locations1[1]-1);
			mapped_contig[tab_locations1[2]-tab_locations1[1]-1]='\0';
			mapping_loc=strtol(tab_locations1[2]+1,NULL,10);
			target_size=strtol(tab_locations2[2]+1,NULL,10);
			if(mapping_loc>target_size) //1st read is read 2 of pair (read order flipped in pair during fastq generation)
			{
				mapped_orientation='-';
				mapping_loc2=mapping_loc;
				mapping_loc=target_size;
				target_size=strtol(tab_locations2[7]+1,NULL,10)+2; //mrFAST insert size = mapping_loc_3'read - mapping_loc_5'read + read_length - 2
				strncpy(cigar,tab_locations2[4]+1,(tab_locations2[5]-tab_locations2[4]-1));
				cigar[tab_locations2[5]-tab_locations2[4]-1]='\0';
				strncpy(cigar2,tab_locations1[4]+1,(tab_locations1[5]-tab_locations1[4]-1));
				cigar2[tab_locations1[5]-tab_locations1[4]-1]='\0';
				strncpy(read1seq,tab_locations2[8]+1,76);
        strncpy(read2seq,tab_locations1[8]+1,76);
				strncpy(original_quality,tab_locations2[9]+1,76);
      	strncpy(original_quality2,tab_locations1[9]+1,76);
				strncpy(md,strstr(line2,"MD")+strlen("MD")+3,300);
				strncpy(md2,strstr(line,"MD")+strlen("MD")+3,300);
			}
			else //1st read is read 1 of pair (read order flipped in pair during fastq generation)
			{
				mapped_orientation='+';
				mapping_loc2=target_size;
				target_size=strtol(tab_locations1[7]+1,NULL,10)+2;
				strncpy(cigar,tab_locations1[4]+1,(tab_locations1[5]-tab_locations1[4]-1));
				cigar[tab_locations1[5]-tab_locations1[4]-1]='\0';
    		strncpy(cigar2,tab_locations2[4]+1,(tab_locations2[5]-tab_locations2[4]-1));
				cigar2[tab_locations2[5]-tab_locations2[4]-1]='\0';
				strncpy(read1seq,tab_locations1[8]+1,76);
        strncpy(read2seq,tab_locations2[8]+1,76);
				strncpy(original_quality,tab_locations1[9]+1,76);
      	strncpy(original_quality2,tab_locations2[9]+1,76);
				strncpy(md,strstr(line,"MD")+strlen("MD")+3,300);
      	strncpy(md2,strstr(line2,"MD")+strlen("MD")+3,300);
			}

			//determine if mapping location corresponds to a valid MIP target
			int wiggle_room=4; //wiggle_room = number of bases actual mapping location can differ from mip target location for a read to be considered to be mapping to a mip target
			int mapped_mip=-1,contig_num;
			for(j=0;j<num_mip_targets;j++)
			{
				contig_num=-1;
				for(k=0;k<max_num_contigs;k++)
				{
					if(strncmp(mips[j].contig_names[k],mapped_contig,17)==0)
					{
						contig_num=k;
						break;
					}
				}
				if((contig_num>-1)&&((mapping_loc-wiggle_room)<=mips[j].contig_start_coords[contig_num])&&((mapping_loc+wiggle_room)>=mips[j].contig_start_coords[contig_num]))
				{
					if(mips[j].orientation==mapped_orientation)
					{
						mapped_mip=j;
						break;
					}
				}
			}
			if(mapped_mip==-1) //mapping location does not correspond to a valid MIP target
			{
				continue;
			}

			//determine if mapped mip target sequence has internal indels between paralogs and ensure insert size is consistent with the expectation for MIPs (152 bp)
			int strict_insert_wiggle_room=2,indel_insert_wiggle_room=10; //insert_wiggle_room = number of bp mapped insert size can differ from expected insert size
			int sum=0;
			int has_indel=0;
			for(k=0;k<152;k++)
			{
				if(mips[mapped_mip].important_bases[k]>0)
				{
					sum++;
				}
			}
			if(sum>=20) 
			//mip target sequence clearly has indel; 20+ paralog-variant bases are due to the fact that aligned sequences were not used to determine these
			//assumption of no indels may be violated for exon-targeting or single-SUN targeting MIPs which were designed blind to the alignment between paralogs
			{
				has_indel=1;
				if((target_size<(152-indel_insert_wiggle_room))||(target_size>(152+indel_insert_wiggle_room)))
				{
					continue; //insert size of mapped reads is outside of expected range for a MIP target 
				}
			}
			else //low number of paralog-variant bases suggests an indel within mip target sequence is not likely
			{
				if((target_size<(152-strict_insert_wiggle_room))||(target_size>(152+strict_insert_wiggle_room)))
      	{
        	continue; //insert size of mapped reads is outside of expected range for a MIP target
      	}
			}

			//initialize array of quality values to store quality scores of read bases corresponding to reference bases 1-152 of MIP target sequence
			//also initialize array of 0s and 1s specifying if a base mapped to a particular reference position is mismatched and a sequence array
			for(i=0;i<152;i++)
			{
				sequence[i]='N';
				quality[i]='!'; //these correspond to contig reference bases 1-152, ! is lower than the lowest actual quality possible, #
				is_base_mm[i]=0;
			}
			//parse CIGAR string and MD tag for reads with function calls
			int inconsistent_reads=0;
			inconsistent_reads+=parse_cigar_and_md(cigar,md,read1seq,original_quality,mapping_loc,mips[mapped_mip].contig_start_coords[contig_num],sequence,quality,is_base_mm);
			inconsistent_reads+=parse_cigar_and_md(cigar2,md2,read2seq,original_quality2,mapping_loc2,mips[mapped_mip].contig_start_coords[contig_num],sequence,quality,is_base_mm);
			
			if(inconsistent_reads)
			{
				continue;	//forward and reverse reads had different base calls for the same base
			}

			int passed_quality=1;
			for(k=0;k<152;k++) //look at quality scores of paralog-variant bases to ensure high quality scores (> Phred 30) at these bases
			{
				if(mips[mapped_mip].important_bases[k]>0)
				{
					if(quality[k]<'?')
					{
						passed_quality=0;
					}
				}
			}
			sum=0;
			if(has_indel) //if there are indel differences between paralogs, ensure the read has an average base quality > Phred 30
			{
				passed_quality=1;
				for(k=0;k<152;k++)
				{
					sum+=quality[k];
				}
				if((sum/k)<'?')
				{
					passed_quality=0;
				}
			}

			if(!(passed_quality))
			{
				continue; //quality score at a paralog-variant base was < Phred 30
			}
			int mismatch_at_sun=0;
			for(k=0;k<152;k++)
			{
				if((mips[mapped_mip].important_bases[k]*is_base_mm[k])>1)
				{
					mismatch_at_sun=1;
				}
			}
			if(mismatch_at_sun)
			{
				//printf("Possible mismapping (mismatch or deletion at SUN position) for read pair:\n%s\n%s\n",line,line2);
				continue; //mismatch or deletion at SUN position indicates a possible mismapping
			}

			//determine if paralog mapped to can be distinguished (has specificity) or not; if not, change the contig mapped to the "nonspecific contig"
			if(mips[mapped_mip].specificity[contig_num]!=1) //mapping lacks specificity
      {
        contig_num=max_num_contigs;
      }
				
			//determine if molecular tag has been seen before at the MIP target mapped to; if never seen before, add it to the tag information data as a new entry
			int new_tag=1;
			int storage_full=0;
			for(y=0;y<5000;y++)
			{
				if(data->tags[mapped_mip][contig_num][y][0]=='N')
				{
					break; //if this occurs it means you have gone through each distinct tag sequence already entered for this MIP-paralog combination without finding a match; thus, your tag sequence is novel
				}
				if((strncmp(tagseq,data->tags[mapped_mip][contig_num][y],taglength)==0))
				{
					new_tag=0; //match found; tag sequence is not novel for this MIP-paralog combination
          break;
				}
				if(y==4999)
				{
					storage_full=1;
				}	
			}
			if((new_tag)&&(!(storage_full))) //tag sequence has not yet been observed and storage is not yet full
      {   
				strncpy(data->tags[mapped_mip][contig_num][y],tagseq,taglength);
      }

			//increment appropriate count if molecular tag has not been seen before at that MIP target
			if(new_tag)
			{
				(data->counts[mapped_mip][contig_num])++;
			}
		}
		gzclose(in);
	}

	//print counts to standard output
	char individual[51];
	strncpy(individual,input_file_name,(strstr(input_file_name,"_FS")-input_file_name));
	individual[(strstr(input_file_name,"_FS")-input_file_name)]='\0';
	long indiv=strtol(*(argv+4),NULL,10);
	if(indiv==1) //print header before data for first individual but not before data for other individuals
	{
		printf("Individual\tTarget_Sequence\tTarget_Coordinate\tMip_Type\t");
		for(k=0;k<max_num_contigs;k++)
		{
			printf("Paralog_%d_Count\t",k+1);
		}
		printf("Nonspecific_Count\n");
	}
	for(j=0;j<num_mip_targets;j++)
	{
		printf("%s\t%s\t%ld\t%c\t",individual,mips[j].master_target_name,mips[j].master_coordinate,mips[j].mip_type);
		for(k=0;k<max_num_contigs;k++)
		{
			printf("%ld\t",data->counts[j][k]);
		}
		printf("%ld\n",data->counts[j][max_num_contigs]);
	}
	
	//clean up and exit
	free(mips);
	free(data);
	return 0;
}

//function to parse through CIGAR string and MD tag corresponding to a single read
//inputs: pointer to cigar string, pointer to md string, pointer to array of read quality scores, mapping location, expected mapping location, pointer to array of final (aligned to contig reference) quality scores,
//pointer to array of integers indicating whether a particular base is mismatched
//output: no output, but the values of the final quality array and mismatched bases array are altered appropriately

//the algorithm to parse CIGAR string is fairly straightforward:
//	first array index = mapping location - theoretical mapping location
//	for mapped bases (M), set array value to corresponding mapped base quality and increment base quality and array index
//	for inserted bases relative to reference, increment base quality (skip over base quality, moving to value at next read base) without changing array index
//	for deleted bases relative to reference, increment array index without changing base quality position
//the algorithm to parse MD tag involves getting the number of matches, getting the next character in the string (either A, C, G, T, or ^), and performing
//the next action accordingly; it has been shown to work on a wide spectrum of test cases
//the CIGAR string and MD tag must be parsed simultaneously because they inform each other's interpretation, see http://davetang.org/muse/2011/01/28/perl-and-sam/

int parse_cigar_and_md(char*cigar_string,char*md_string,char*read_sequence,char*original_quality_array,long mapping_location,long expected_mapping_location,char*sequence_array,char*quality_array,int*is_base_mm_array)
{
	long index=mapping_location-expected_mapping_location; //starting index of final quality array and mismatch array
	int read_index=0; //starting index of the read quality array
	int inconsistency=0;
	long num_bases,num_matches;
	int md_offset=0;
	char alignment_type,ch;
  char dummy[4];
	int i=0,j,k;
	int from_ins=0;
	while(cigar_string[i]!='\0') //while some CIGAR string left to process
	{
		num_bases=strtol(cigar_string+i,NULL,10); //calculate number of bases mapped (M), inserted (I), or deleted (D)
		while(!(isalpha(cigar_string[i])))
    {
    	i++;
    }
    alignment_type=cigar_string[i]; //calculate whether bases mapped (M), are inserted (I), or are deleted (D)
		i++; //increment i so next part of cigar string will be either a number or the null character
		if(alignment_type=='M') //read bases mapped
		{
			if(!from_ins)
			{
				k=0;
				num_matches=strtol(md_string+md_offset,NULL,10); //determine how many mapped bases are matches from the MD string
				sprintf(dummy,"%ld\0",num_matches);
      	md_offset+=strlen(dummy);
      	ch=md_string[md_offset]; //using MD string, determine whether first base after last matching base is a substitution (A,C,G,T) or 'deleted' (^) compared to the reference
			}
			from_ins=0;
			for(j=0;j<num_bases;j++,k++)
			{
				if(k==num_matches) //at this point index corresponds to the array location after the final matching base
				{
					if(isalpha(ch))
					{
						md_offset++;
						if((index>=0)&&(index<152))
						{
							is_base_mm_array[index]=1;
							if((sequence_array[index]!='N')&&(sequence_array[index]!=read_sequence[read_index])) //if analyzing 2nd read and 1st read base does not match 2nd read base
							{
								inconsistency=1;
							}
							sequence_array[index]=read_sequence[read_index];
						}
						num_matches=strtol(md_string+md_offset,NULL,10);
      			sprintf(dummy,"%ld\0",num_matches);
      			md_offset+=strlen(dummy);
      			ch=md_string[md_offset];
						k=-1; //since k is incremented at the end of the loop and we want new k=0
					}	
				}
				if((quality_array[index]=='!')&&(index>=0)&&(index<152))
				{
					quality_array[index]=original_quality_array[read_index];
					sequence_array[index]=read_sequence[read_index];
				}
				else if((index>=0)&&(index<152)) //there is some overlap between 1st read and 2nd read due to deletion(s); take the highest quality if there is agreement on the base call
				{
					if(sequence_array[index]==read_sequence[read_index])
					{
						quality_array[index]=(quality_array[index]>original_quality_array[read_index])?quality_array[index]:original_quality_array[read_index];
					}
					else
					{
						inconsistency=1;
					}
				}
				read_index++; //increment read index
				index++;	//increment final index
			}
		}
		if(alignment_type=='I') //read bases are inserted compared to the contig reference
		{
			from_ins=1;
			for(j=0;j<num_bases;j++)
			{
				read_index++; //increment read index only
			}
		}
		if(alignment_type=='D') //read lacks bases compared to the contig reference (some reference bases are deleted)
		{
			from_ins=0;
			md_offset+=(num_bases+1); //the '+1' is for the '^' character
			for(j=0;j<num_bases;j++)
			{
				if((quality_array[index]=='!')&&(index>=0)&&(index<152))
				{
					is_base_mm_array[index]=1;
					sequence_array[index]='-';
					if((read_index>=1)&&(read_index<76))
					{
						quality_array[index]=(original_quality_array[read_index-1]+original_quality_array[read_index])/2; //quality values for 'deleted' base positions should be the average of the flanking quality values in the read
					}
				}
				else if((index>=0)&&(index<152)&&(read_index>=1)&&(read_index<76))
				{
					quality_array[index]=(quality_array[index]>((original_quality_array[read_index-1]+original_quality_array[read_index])/2))?quality_array[index]:(original_quality_array[read_index-1]+original_quality_array[read_index])/2;
				}
				index++; //increment final index only
			}
		}		
	}
	return inconsistency;
}
