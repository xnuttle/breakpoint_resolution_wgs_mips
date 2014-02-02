//Xander Nuttle
//qseq_to_fastq_for_mrfast_tagged_split.c
//Takes three gzipped qseq text files containing first reads, index reads, and second reads, respectively, and generates gzipped fastq files as output.
//Only adds reads having a barcode sequence perfectly matching a known barcode reverse complement sequence to the output fastq files.
//Generates 1 or more gzipped fastq files per individual-barcode combination, demultiplexing the reads in the process.
//Assumes MIPs have been molecularly tagged at one location just internal to the extension arm to allow for counting of individual capture events.
//
//Call: ./qseq_to_fastq_for_mrfast_tagged_split reads1_qseq.txt.gz reads2_qseq.txt.gz reads3_qseq.txt.gz (int)trimmed_read_length (int)molecular_tag_length (int)file_set_number (long)max_num_reads_per_output_file barcode_key_file
//
//This program converts gzipped qseq text files into gzipped fastq files, which serve as input to mappers such as BWA, mrFAST, and mrsFAST.
//Qseq text files containing paired-end sequence data are naturally grouped in threes - the first file contains the first sequence reads from
//read pairs, the second file contains index reads that associate read pairs with a particular sample, and the third file contains the second
//sequence reads from read pairs. Sequence data for each read pair is stored on the same line in each of these files. For example, line 1 of 
//the first file contains data for the first read of the first read pair, line 1 of the second file contains data for the index read corresponding
//to the first read pair, and line 1 of the third file contains data for the second read of the first read pair.
//
//This program performs several functions:
//	-converts a set of three gzipped qseq text files into multiple gzipped fastq files, splitting the reads by sample based on barcode sequences,
//	 with each read in the gzipped fastq files tagged with the molecular tag
//	 (each output gzipped fastq output file will contain a maximum of a user-specified number of reads [500,000 or fewer is optimal for mrfast input files])
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
  barcodekey=fopen(*(argv+8),"r");
  fgetpos(barcodekey,&start_bcfile);
  fscanf(barcodekey,"%s %s",dummystr[0],dummystr[1]);
  fsetpos(barcodekey,&start_bcfile);
  index_read_length=strlen(dummystr[1]); //barcode length
	char sample_names[384][51];
	char barcodes[384][index_read_length+1]; //we have 384 barcodes, so the maximum possible number of individuals per experiment is 384
	int num_indivs;
	int indiv=0;
	while(fscanf(barcodekey,"%s %s",sample_names[indiv],barcodes[indiv])==2)
	{
		indiv++;
	}
	num_indivs=indiv;
	fclose(barcodekey);

	//setup input files
	gzFile*in1,*in2,*in3;

	//read in file set number from the command line
	long fileset; //specifying the number of the file set the program is working on allows parallelization on the cluster
	fileset=strtol(*(argv+6),NULL,10);
	
	//setup output files
	char outname[51];
	gzFile*outfiles[num_indivs];
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		sprintf(outname,"%s_FS%ld_1.fastq.gz\0",sample_names[indiv],fileset);
		outfiles[indiv]=gzopen(outname,"w");
	}

	//setup variables: specify read length, index read length, and desired number of sequence reads per fastq file
	long trimmed_read_length; //76 suggested for MIPs (76 bp is the minimal value to cover all targeted bases (112 bp) + both hybridization arms (20 bp each)), 100 bp recommended for indel calling
	long tag_length;
	long reads_per_fastq=strtol(*(argv+7),NULL,10); //500,000 is the recommended number of total reads per fastq file for optimal performance of mrfast and mrsfast mappers
	int output_file_nums[num_indivs];
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		output_file_nums[indiv]=1;
	}

	//read in trimmed read length value from the command line
	trimmed_read_length=strtol(*(argv+4),NULL,10);

	//read in molecular tag length value from the command line
	tag_length=strtol(*(argv+5),NULL,10);
	
	//setup variables
	in1=gzopen(*(argv+1),"r");
  in2=gzopen(*(argv+2),"r");
  in3=gzopen(*(argv+3),"r");
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
	char tag_sequence[tag_length+1];
  tag_sequence[tag_length]='\0';
	
	//read in gzipped qseq text files line by line, trim sequences, and print output
	while((line[0]=gzgetc(in3))!=-1) //a -1 return value for gzgetc() is equivalent to EOF
	{		
		for(indiv=0;indiv<num_indivs;indiv++)
		{
			if(reads_output[indiv]>=reads_per_fastq)
			{
				gzclose(outfiles[indiv]);
				output_file_nums[indiv]++;
				sprintf(outname,"%s_FS%ld_%d.fastq.gz\0",sample_names[indiv],fileset,output_file_nums[indiv]);
				outfiles[indiv]=gzopen(outname,"w");
				reads_output[indiv]=0;
			}
		}
		i=1;
		m=0;
  	
		//process data for second sequence
		while((line[i]=gzgetc(in3))!='\n')
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
		strncpy(sequence,tab_locations[7]+1+tag_length,trimmed_read_length);	
		strncpy(quality,tab_locations[8]+1+tag_length,trimmed_read_length);
		strncpy(tag_sequence,tab_locations[7]+1,tag_length);
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
		while((line[i]=gzgetc(in2))!='\n')
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

		//output second sequence to gzipped fastq file
		if((indiv!=-1)&&(!(strchr(tag_sequence,'N')))&&(strstr(tag_sequence,"AAAAA")==NULL)&&(strstr(tag_sequence,"CCCCC")==NULL)&&(strstr(tag_sequence,"GGGGG")==NULL)&&(strstr(tag_sequence,"TTTTT")==NULL))
		{
			//barcode read perfectly matches a known barcode, molecular tag sequence does not contain any ambiguous bases, and molecular tag does not contain any homopolymer > 4 nt
			gzprintf(outfiles[indiv],"@%s:%ld:%ld:%ld:%ld:%c$%s/1\n%s\n+\n%s\n",machine_name,flow_cell_lane,tile_number,cluster_x_coord,cluster_y_coord,passed_chastity_filter,tag_sequence,sequence,quality);
		}
	
		//process data for first sequence
		i=0;
		m=0;
		while((line[i]=gzgetc(in1))!='\n')
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

		//output first sequence to gzipped fastq file
    if((indiv!=-1)&&(!(strchr(tag_sequence,'N')))&&(strstr(tag_sequence,"AAAAA")==NULL)&&(strstr(tag_sequence,"CCCCC")==NULL)&&(strstr(tag_sequence,"GGGGG")==NULL)&&(strstr(tag_sequence,"TTTTT")==NULL))
		{
			gzprintf(outfiles[indiv],"@%s:%ld:%ld:%ld:%ld:%c$%s/2\n%s\n+\n%s\n",machine_name,flow_cell_lane,tile_number,cluster_x_coord,cluster_y_coord,passed_chastity_filter,tag_sequence,sequence,quality);
			reads_output[indiv]+=2;
		}
	}
	
	//clean up and exit
	gzclose(in1);
  gzclose(in2);
  gzclose(in3);
	for(indiv=0;indiv<num_indivs;indiv++)
	{
		gzclose(outfiles[indiv]);
	}
	return 0;
}
