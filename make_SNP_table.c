//Xander Nuttle
//make_SNP_table.c
//Call: make_SNP_table prox_aligned.fasta dist_aligned.fasta <other_aligned.fasta>
//
//This program generates a SNP table for MIP design over paralogous sequences. 'SNPs' are defined in this case as positions of the master sequence where the alignment
//includes a mismatch or a gap, as well as positions of the master sequence adjacent to an alignment gap in the master sequence (as gap positions are not truly
//positions in the master sequence). The SNP table output by this program will inform the MIP design program to flag MIPs with arms overlapping these alignment
//positions as having SNPs. Thus, such MIPs, which would likely not hybridize equally well to all paralogous sequences of interest, can be eliminated from further
//consideration for design. Fasta files with '-'s containing aligned sequences serve as inputs to this program. These files should be provided as command line
//arguments, with the "master" sequence the first file input to this program.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char*argv[])
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
  sprintf(outfilename,"%s.snptable",mastername);
  FILE*out=fopen(outfilename,"w");
	
	//compare characters in aligned paralogous sequences and output any positions in master sequence that are not identical between all paralogs
	long delbases=0;
	long k;
	int n,diff;	
	for(k=0;k<j;k++)
	{
		if(seqs[0][k]=='-')
		{
			delbases++;
			if(diff<1)
			{
				fprintf(out,"585 chr%s_shared %ld %ld SNP 0 + A A A/C genomic unknown unknown 0 0 unknown exact 1\n",mastername,k-delbases,k+1-delbases);
				diff=1;
			}
		}
		else
		{
			diff=0;
			for(i=0;((i<(argc-2))&&(diff<1));i++)
			{
				for(n=1;((n<(argc-1))&&(diff<1));n++)
				{
					if(seqs[i][k]!=seqs[n][k])
					{
						fprintf(out,"585 chr%s_shared %ld %ld SNP 0 + A A A/C genomic unknown unknown 0 0 unknown exact 1\n",mastername,k-delbases,k+1-delbases);
						diff=1;
					}
				}
			}
		}
	}

	//cleanup and exit
	for(i=0;i<(argc-1);i++)
  {
    fclose(files[i]);
    free(seqs[i]);
  }
  fclose(out);
  return 0;
}
