//Xander Nuttle
//brkpt_mipmaker.c
//Call: brkpt_mipmaker prox_aligned.fasta dist_aligned.fasta <other_aligned.fasta>
//
//This program identifies all regions from a multiple alignment of paralogous sequences where MIP design should be attempted. Such regions are defined by a 112 bp
//region of the alignment where at least one aligned sequence is distinct and where there are no gaps in the alignment, flanked on both sides by at least 16 bp of
//sequence that is identical between all sequences in the alignment. The multiple alignment is input via fasta files with '-'s containing the aligned sequences.
//These files should be provided as command line arguments, with the "master" sequence the first file input to this program.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int arms_identical(int numseqs,char*seqs[numseqs],long start,long armlength);
int target_uniq(int numseqs,char*seqs[numseqs],long start,long targlength);

int main(int argc,char*argv[])
{
	//open input files and read and store master sequence name
	FILE*files[argc-1];
  int i;
	long j;
  for(i=0;i<(argc-1);i++)
  {
    files[i]=fopen(*(argv+1+i),"r");
  }
	char mastername[51];
	for(i=0;i<51;i++)
	{
		mastername[i]='\0';
	}
	strncpy(mastername,*(argv+1),strrchr(*(argv+1),'_')-(*(argv+1)));

	//determine length of master sequence, including '-'s in the alignnment
	long alignsize=0;
	char ch;
	while((ch=getc(files[0]))!='\n')
  	continue;
	while((ch=getc(files[0]))!=EOF)
  {
  	if((isalpha(ch))||(ch=='-'))
    	alignsize++;
	}
	fclose(files[0]);
	files[0]=fopen(*(argv+1),"r");

	//setup arrays to store aligned sequences
	char*seqs[(argc-1)];
	for(i=0;i<(argc-1);i++)
	{
		seqs[i]=(char*)malloc((alignsize+1)*sizeof(char));
		if(seqs[i]==NULL)
		{
			printf("Memory allocation failed!\n");
    	return 1;
		}
	}

	//read sequences into character arrays
	char base;
	for(i=0;i<(argc-1);i++)
	{
		while((base=getc(files[i]))!='\n')
			continue;
		j=0;
		while((base=getc(files[i]))!=EOF)
		{
			if(!isspace(base))
    	{
      	seqs[i][j]=base;
      	j++;
    	}
		}
		seqs[i][j]='\0';
	}
	
	//setup output file
	char outfilename[76];
	outfilename[75]='\0';
	sprintf(outfilename,"%s_regions_for_design.bed",mastername);
	FILE*out=fopen(outfilename,"w");

	//scan character arrays for suitable target regions
	long target_size=112,arm_length=16;
	long delbases=0;
	int has_gap=0;
	char*targseqs[(argc-1)];
  for(i=0;i<(argc-1);i++)
  {
    targseqs[i]=(char*)malloc((target_size+1)*sizeof(char));
    if(seqs[i]==NULL)
    {
      printf("Memory allocation failed!\n");
      return 1;
    }
  }
	for(j=0;j<=strlen(seqs[0])-(target_size+(2*arm_length));j++)
	{
		if(seqs[0][j]=='-')
		{
			delbases++;
		}
		if(arms_identical((argc-1),seqs,j,arm_length))
    {
			if(arms_identical((argc-1),seqs,(j+arm_length+target_size),arm_length))
      {
				if(target_uniq((argc-1),seqs,(j+arm_length),target_size))
				{
					has_gap=0;
					for(i=0;i<(argc-1);i++)
					{
						strncpy(targseqs[i],seqs[i]+j+arm_length,target_size);
						targseqs[i][target_size]='\0';
						if(strchr(targseqs[i],'-')!=NULL)
						{
							has_gap=1;
						}
					}
					if(!(has_gap))
					{
						fprintf(out,"chr%s_shared\t%ld\t%ld\n",mastername,j+arm_length-delbases,j+arm_length+target_size-delbases);
					}
				}
			}
		}
	}

	//clean up and exit
	for(i=0;i<(argc-1);i++)
	{
		fclose(files[i]);
		free(seqs[i]);
		free(targseqs[i]);
	}
	fclose(out);
	return 0;
}

int arms_identical(int numseqs,char*seqs[numseqs],long start,long armlength)
{
  int identical=1;
  int seq1,seq2;
  for(seq1=0;seq1<(numseqs-1);seq1++)
  {
    for(seq2=(seq1+1);seq2<(numseqs);seq2++)
    {
			if(strncmp(seqs[seq1]+start,seqs[seq2]+start,armlength)!=0)
      {
				identical=0;
      }
    }
  }
  return identical;
}

int target_uniq(int numseqs,char*seqs[numseqs],long start,long targlength)
{
  int uniq=0;
  int seq1,seq2;
  for(seq1=0;seq1<(numseqs-1);seq1++)
  {
    for(seq2=(seq1+1);seq2<(numseqs);seq2++)
    {
      if(strncmp(seqs[seq1]+start,seqs[seq2]+start,targlength)!=0)
      {
        uniq=1;
      }
    }
  }
  return uniq;
}
