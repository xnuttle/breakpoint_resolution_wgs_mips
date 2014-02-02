//Xander Nuttle
//get_contig_SUN_coords.c
//Call: get_contig_SUN_coords prox_aligned.fasta dist_aligned.fasta <other_aligned.fasta> alignment_base1_is_contig_base_N.coordconvert
//
//SUN characterization is important for identifying single-nucleotide markers that can be leveraged to refine breakpoint locations within segmental duplications.
//This program identifies SUNs from fasta files containing '-'s output from Jalview or another alignment editor and outputs their coordinates, both with respect
//to the alignment (specifically, with respect to the aligned sequence in the first file input to this program) and with respect to contig fasta files lacking '-'s.
//All output coordinates are in base 1. This program can accomodate sequences paralogous to breakpoint-associated segmental duplications. These sequences should
//be aligned together with the breakpoint-associated segmental duplications, and fasta files with '-'s containing their aligned sequences should be included as
//additional command line inputs to this program. The final command line input should be a tab-delimited three column text file, with the first column having the
//names of the contigs corresponding to each aligned sequence, the second column having the positions (in base 1 coordinates) of the bases in those contigs
//corresponding to the first base in each aligned sequence, and the third column having the orientations of the aligned sequences to their corresponding contigs.
//
//The first column of the output file is the name of the contig having the corresponding SUN, the second column is the base 1 coordinate of the SUN with respect
//to its corresponding contig sequence, and the third column is the base 1 coordinate of the SUN with respect to the "master" sequence, the aligned sequence in the
//first file input to this program. The fourth and fifth columns contain placeholders previously used in MIP design retained here only because a later program expects
//them to be present.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc,char*argv[])
{
	//open all input files, read and store contig names, setup array to keep track of deleted bases in each aligned sequence
	FILE*files[argc-1];
	int i,j;
	for(i=0;i<(argc-1);i++)
	{
		files[i]=fopen(*(argv+1+i),"r");
	}
	char contignames[i-1][51];
	char orientation[i-1];
	long alignstarts[i-1],aligncoords[i-1],delbases[i-1];
	for(i=0;i<(argc-2);i++)
	{
		for(j=0;j<51;j++)
		{
			contignames[i][j]='\0';
		}
		delbases[i]=0;
		fscanf(files[argc-2],"%s %ld %c",contignames[i],&(alignstarts[i]),&(orientation[i]));
	}

	//setup output file
	FILE*out;
	out=fopen("contig.suns","w");

	//go through alignment output fasta files and output information about all SUNs
	char base[i+1];
  base[i]='\0';
	for(i=0;i<(argc-2);i++)
	{
		while((base[i]=getc(files[i]))!='\n')
			continue;
	}
	for(i=0;i<(argc-2);i++)
  {
    aligncoords[i]=1;
  }
	while((base[0]=getc(files[0]))!=EOF)
	{
		for(i=1;i<(argc-2);i++)
		{
			base[i]=getc(files[i]);
		}
		for(i=0;i<(argc-2);i++)
		{
			while(isspace(base[i]))
			{
      	base[i]=getc(files[i]);
    	}
		}
		for(i=0;i<(argc-2);i++)
    {
      if(base[i]=='-')
			{
				delbases[i]++;
			}
    }
		for(i=0;i<(argc-2);i++)
    {
      if((strchr(base,base[i])==strrchr(base,base[i]))&&(strpbrk(base,"YRWSKMDVHBN-")==NULL)) //SUNs
      {
				if(orientation[i]=='+')
				{
					fprintf(out,"%s\t%ld\t%ld\t%lf\t%s\n",contignames[i],aligncoords[i]-delbases[i]+alignstarts[i]-1,aligncoords[0]-delbases[0],10.0,"no_masking_near_SUN");
				}
				else
				{
					fprintf(out,"%s\t%ld\t%ld\t%lf\t%s\n",contignames[i],alignstarts[i]-(aligncoords[i]-delbases[i])+1,aligncoords[0]-delbases[0],10.0,"no_masking_near_SUN");
				}
			}
		}
		for(i=0;i<(argc-2);i++)
    {
        aligncoords[i]++;
    }
	}

	//cleanup and exit
	for(i=0;i<(argc-1);i++)
	{
		fclose(files[i]);
	}
	fclose(out);
	return 0;
}
