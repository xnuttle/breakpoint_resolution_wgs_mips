//Xander Nuttle
//get_sunk_reads.c
//Call: get_sunk_reads sunks_file text_file_with_names_of_gzipped_qseq_files
//
//This program prints to a file all reads having a perfect 36 bp match to a breakpoint-informative SUNK, as well as associated information (in fastq format).
//These reads can then be analyzed to calculate read depth at the SUNs associated with these SUNKs. Visualization of this data should allow breakpoint
//locations to be inferred.

#include<stdio.h>
#include<zlib.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>

void revcomp(char string[]);
unsigned long bsi(char*kmer,char**pointer_array,unsigned long num_kmers);
int compfun(const void*p1,const void*p2);
int sunklength=36;

int main(int argc,char*argv[])
{
	//read in SUNKs, calculate their reverse complements, and sort SUNKs and their reverse complements
	long i,j;
	FILE*sunksfile;
	int numsunks=0,start;
	sunksfile=fopen(*(argv+1),"r");
  fpos_t pos;
	fgetpos(sunksfile,&pos);	
	char seqname[51];
	seqname[50]='\0';
  while(fscanf(sunksfile,"%s %ld %s",seqname,&start,seqname)==3)
  {
			numsunks++;
	}
	fsetpos(sunksfile,&pos);
	char**sunks;
	sunks=(char**)malloc((2*numsunks+2)*sizeof(char*));
	char firststring[]="!";
	char laststring[]="~";
	if(sunks==NULL)
  {
    printf("Memory allocation failed!\n");
    return 1;
  }
	for(i=0;i<(2*numsunks+2);i++)
	{
		sunks[i]=(char*)malloc((sunklength+1)*sizeof(char));
		if(sunks[i]==NULL)
		{
			printf("Memory allocation failed!\n");
    	return 1;
		}
		for(j=0;j<(sunklength+1);j++)
    {
      sunks[i][j]='\0';
    }
	}
	strncpy(sunks[0],firststring,sunklength);
	i=1;
	while(fscanf(sunksfile,"%s %ld %s",seqname,&start,sunks[i])==3)
  {
		strncpy(sunks[i+1],sunks[i],sunklength);
		revcomp(sunks[i+1]);
		i+=2;
  }
	strncpy(sunks[i],laststring,sunklength);
	fclose(sunksfile);
	qsort(sunks,(2*numsunks+2),sizeof(char*),compfun);

	//setup input gzipped fastq files
	FILE*textfile;
  textfile=fopen(*(argv+2),"r");
  gzFile*fastq;
  char file[51];
	file[50]='\0';

	//setup output file
	FILE*out;
	out=fopen("filtered_reads.fastq","w");	

	//process fastq files and print to output file all reads having a kmer matching a SUNK
	int readlength;
  unsigned long loc;
  char line1[2001];
	char line2[2001];
	char line3[2001];
  char line4[2001];
	while(fscanf(textfile,"%s",file)==1)
	{
  	fastq=gzopen(file,"r");
  	gzgets(fastq,line2,2000);
		for(i=0;i<2001;i++)
  	{
    	line2[i]='\0';
  	}
  	gzgets(fastq,line2,2000);
		gzrewind(fastq);
  	readlength=(strlen(line2)-1);
		for(i=0;i<2001;i++)
  	{
    	line1[i]='\0';
    	line2[i]='\0';
    	line3[i]='\0';
    	line4[i]='\0';
  	}
		char kmer[sunklength+1];
  	kmer[sunklength]='\0';
  	while(gzgets(fastq,line1,2000))
  	{
			gzgets(fastq,line2,2000);
			gzgets(fastq,line3,2000);
			gzgets(fastq,line4,2000);
			for(i=0;i<(readlength-sunklength+1);i++)
			{
				strncpy(kmer,(line2+i),sunklength);
				loc=bsi(kmer,sunks,(2*numsunks+2));
				if(loc!=ULONG_MAX)
				{
					fprintf(out,"%s%s%s%s",line1,line2,line3,line4);
					break;
				}
			}
			if(gzeof(fastq))
			{
				break;
			}
   	}
		gzclose(fastq);
	} 
	
	//cleanup and exit
	for(i=0;i<(2*numsunks+2);i++)
	{
		free(sunks[i]);
	}
	free(sunks);
	fclose(out);
	return 0;
}

void revcomp(char string[])
{
  int a=strlen(string);
  int b,c;
  char base;
  char newstring[a+1];
  for(b=0;b<a;b++)
  {
  newstring[b]=string[a-1-b];
  }
  for(c=0;c<a;c++)
  {
    base=newstring[c];
    switch(base)
    {
    case 'A':   base='T';
          break;
    case 'C':   base='G';
          break;
    case 'G': base='C';
          break;
    case 'T': base='A';
          break;
    }
    newstring[c]=base;
  }
  newstring[a]='\0';
  strcpy(string,newstring);
}

int compfun(const void*p1,const void*p2)
{
  char*const*a1=p1;
  char*const*a2=p2;
  char s1[sunklength+1],s2[sunklength+1];
  strncpy(s1,*a1,sunklength);
  strncpy(s2,*a2,sunklength);
  return(strcmp(s1,s2));
}

unsigned long bsi(char*kmer,char**pointer_array,unsigned long num_kmers)
{
  unsigned long start=0;
  unsigned long end=num_kmers;
  unsigned long index;
  int i;
  while(start<end)
  {
    index=(start+end)/2;
    i=strncmp(kmer,pointer_array[index],sunklength);
    if(i==0) {return index;}
    else if(i>0) {start=index+1;}
    else {end=index-1;}
  }
  index=(start+end)/2;
  i=strncmp(kmer,pointer_array[index],sunklength);
  if(i==0) {return index;}
  return ULONG_MAX;
}
