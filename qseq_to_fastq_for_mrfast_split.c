//Xander Nuttle
//qseq_to_fastq_for_mrfast_split.c
//Takes a file listing sets of gzipped qseq text files containing first reads, index reads, and second reads, respectively, and generates gzipped fastq files as output.
//Only adds reads having a barcode sequence perfectly matching a known barcode reverse complement sequence to the output fastq files.
//Generates 1 gzipped fastq file per individual-barcode combination, demultiplexing the reads in the process.
//
//Call: ./qseq_to_fastq_for_mrfast_split text_file_with_names_of_gzipped_qseq_files (int)trimmed_read_length (long)max_num_reads_per_output_file barcode_key_file
//
//This program converts gzipped qseq text files into gzipped fastq files, which serve as input to mappers such as BWA, mrFAST, and mrsFAST.
//Qseq text files containing paired-end sequence data are naturally grouped in threes - the first file contains the first sequence reads from
//read pairs, the second file contains index reads that associate read pairs with a particular sample, and the third file contains the second
//sequence reads from read pairs. Sequence data for each read pair is stored on the same line in each of these files. For example, line 1 of 
//the first file contains data for the first read of the first read pair, line 1 of the second file contains data for the index read corresponding
//to the first read pair, and line 1 of the third file contains data for the second read of the first read pair.
//
//This program performs several functions:
//	-converts gzipped qseq text files listed in a text file in sets of three into multiple gzipped fastq files, splitting the reads by sample based on barcode sequences
//	 (each output gzipped fastq output file will contain a maximum of a user-specified number of reads [500,000 is the optimal number for mrfast input files])
//	-trims reads, retaining only a user-specified number of bases of each sequence read (from the 5' ends of reads, eliminating lower quality
//	 bases at the 3' ends of reads)
//	-converts ambiguous bases in gzipped qseq text files (specified by ".") to "N"
//
//This program does not filter reads based on quality. The assumption is that low-quality reads will not map within the correct MIP insert size range
//(148-156 bp, allowing for 4 inserted or deleted bp) and thus can be eliminated in the downstream analysis. This approach is advantageous because it 
//doesn't result in usable reads being thrown out based on not passing a stringent quality filter. It also prevents only one read from a read pair being
//input to the mapper. Of course, individual base qualities are taken into account in the final paralog-specific read counting analysis.

#include<stdio.h>
#include<zlib.h>
#include<string.h>
#include<stdlib.h>

int main(int argc,char*argv[])
{
	//read in barcode key file and determine barcode length and the number of individuals in the experiment
	FILE*barcodekey;
	char dummystr[2][51];
  fpos_t start_bcfile;
  int index_read_length;
  barcodekey=fopen(*(argv+4),"r");
  fgetpos(barcodekey,&start_bcfile);
  fscanf(barcodekey,"%s %s",dummystr[0],dummystr[1]);
  fsetpos(barcodekey,&start_bcfile);
  index_read_length=strlen(dummystr[1]); //barcode length
	char sample_names[384][51];
	char barcodes[384][index_read_length+1]; //we have 384 barcodes, so the maximum possible number of individuals per experiment is 384; each barcode is 8 bp
	int num_indivs;
	int indiv=0;
	while(fscanf(barcodekey,"%s %s",sample_names[indiv],barcodes[indiv])==2)
	{
		indiv++;
	}
	num_indivs=indiv;
	fclose(barcodekey);

	//setup input files
	FILE*textfile;
	textfile=fopen(*(argv+1),"r");
	gzFile*in[3];
	char input_file_names[3][51];
	
	//setup output files
	char outname[51];
	gzFile*outfiles[num_indivs];
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		sprintf(outname,"%s_1.fastq.gz\0",sample_names[indiv]);
		outfiles[indiv]=gzopen(outname,"w");
	}

	//setup variables: specify read length, index read length, and desired number of sequence reads per fastq file
	long trimmed_read_length; //76 suggested for MIPs (76 bp is the minimal value to cover all targeted bases (112 bp) + both hybridization arms (20 bp each)), 100 bp recommended for indel calling
	int read_length=151;
	long reads_per_fastq=strtol(*(argv+3),NULL,10); //500,000 is the recommended number of total reads per fastq file for optimal performance of mrfast and mrsfast mappers
	int output_file_nums[num_indivs];
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		output_file_nums[indiv]=1;
	}

	//read in trimmed read length value from the command line
	trimmed_read_length=strtol(*(argv+2),NULL,10);
	
	//setup variables
	char line[501];
	line[500]='\0';
	char machine_name[31];
	machine_name[30]='\0';
	long flow_cell_lane,tile_number,cluster_x_coord,cluster_y_coord,pass_chastity_filter;
	long reads_output[num_indivs];
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		reads_output[indiv]=0;
	}
	char passed_chastity_filter;
	int i,m,k;
	char*tab_locations[10];
	char sequence[trimmed_read_length+1],quality[trimmed_read_length+1];
	sequence[trimmed_read_length]='\0';
	quality[trimmed_read_length]='\0';
	char index_sequence[index_read_length+1],index_quality[index_read_length+1];
	index_sequence[index_read_length]='\0';
	index_quality[index_read_length]='\0';
	
	//read in gzipped qseq text files line by line, trim sequences, and print output
	while(fscanf(textfile,"%s%s%s",&(input_file_names[0]),&(input_file_names[1]),&(input_file_names[2]))==3)
	{
		for(m=0;m<3;m++)
		{
			in[m]=gzopen(input_file_names[m],"r");
		}	
		while((line[0]=gzgetc(in[0]))!=-1) //a -1 return value for gzgetc() is equivalent to EOF
  	{
			for(indiv=0;indiv<num_indivs;indiv++)
			{
				if(reads_output[indiv]>=reads_per_fastq)
				{
					gzclose(outfiles[indiv]);
					output_file_nums[indiv]++;
					sprintf(outname,"%s_%d.fastq.gz\0",sample_names[indiv],output_file_nums[indiv]);
					outfiles[indiv]=gzopen(outname,"w");
					reads_output[indiv]=0;
				}
			}
			i=1;
			m=0;
  	
			//process data for first sequence
			while((line[i]=gzgetc(in[0]))!='\n')
  		{
				if(line[i]=='\t')
				{
					tab_locations[m]=line+i;
					m++;
				}
				if(line[i]=='.')
				{
					line[i]='N';
				}
  			i++;
  		}
			strncpy(machine_name,line,tab_locations[0]-line);
			machine_name[tab_locations[0]-line]='\0';
			flow_cell_lane=strtol(tab_locations[1],NULL,10);
			tile_number=strtol(tab_locations[2],NULL,10);
			cluster_x_coord=strtol(tab_locations[3],NULL,10);
			cluster_y_coord=strtol(tab_locations[4],NULL,10);
			strncpy(sequence,tab_locations[7]+1,trimmed_read_length);	
			strncpy(quality,tab_locations[8]+1,trimmed_read_length);
			pass_chastity_filter=strtol(tab_locations[9],NULL,10);		
			if(pass_chastity_filter)
			{
				passed_chastity_filter='Y';
			}
			else
			{
				passed_chastity_filter='N';
			}
		
			//process data for index sequence
			i=0;
			m=0;
			while((line[i]=gzgetc(in[1]))!='\n')
    	{
      	if(line[i]=='\t')
      	{
        	tab_locations[m]=line+i;
        	m++;
      	}
      	if(line[i]=='.')
      	{
        	line[i]='N';
      	}
      	i++;
    	}
			strncpy(index_sequence,tab_locations[7]+1,index_read_length);
			//index quality will not later explicity be taken into account, but it will in the sense that only index sequences with perfect matches to known barcodes will be used

			//ensure barcode sequence perfectly matches a known barcode reverse complement sequence, and if so, identify which individual the read pair corresponds to
			indiv=-1;
			for(k=0;k<num_indivs;k++)
			{
				if(strncmp(index_sequence,barcodes[k],index_read_length)==0)
				{
					indiv=k;
					break;
				}
			}

			//output first sequence to gzipped fastq file
			if(indiv!=-1) //barcode read did not perfectly match any known barcode
			{
				gzprintf(outfiles[indiv],"@%s:%ld:%ld:%ld:%ld:%c/1\n%s\n+\n%s\n",machine_name,flow_cell_lane,tile_number,cluster_x_coord,cluster_y_coord,passed_chastity_filter,sequence,quality);
			}
	
			//process data for second sequence
			i=0;
			m=0;
			while((line[i]=gzgetc(in[2]))!='\n')
    	{
      	if(line[i]=='\t')
      	{
        	tab_locations[m]=line+i;
        	m++;
      	}
      	if(line[i]=='.')
      	{
        	line[i]='N';
      	}
      	i++;
    	}
    	strncpy(machine_name,line,tab_locations[0]-line);
			machine_name[tab_locations[0]-line]='\0';
    	flow_cell_lane=strtol(tab_locations[1],NULL,10);
    	tile_number=strtol(tab_locations[2],NULL,10);
    	cluster_x_coord=strtol(tab_locations[3],NULL,10);
    	cluster_y_coord=strtol(tab_locations[4],NULL,10);
    	strncpy(sequence,tab_locations[7]+1,trimmed_read_length);
    	strncpy(quality,tab_locations[8]+1,trimmed_read_length);
    	pass_chastity_filter=strtol(tab_locations[9],NULL,10);
			if(pass_chastity_filter)
    	{
      	passed_chastity_filter='Y';
    	}
    	else
    	{
      	passed_chastity_filter='N';
    	}

			//output second sequence to gzipped fastq file
    	if(indiv!=-1)
			{
				gzprintf(outfiles[indiv],"@%s:%ld:%ld:%ld:%ld:%c/2\n%s\n+\n%s\n",machine_name,flow_cell_lane,tile_number,cluster_x_coord,cluster_y_coord,passed_chastity_filter,sequence,quality);
				reads_output[indiv]+=2;
			}
		}
		for(m=0;m<3;m++)
		{
			gzclose(in[m]);
		}
	}
	
	//clean up and exit
	fclose(textfile);
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		gzclose(outfiles[indiv]);
	}
	return 0;
}
