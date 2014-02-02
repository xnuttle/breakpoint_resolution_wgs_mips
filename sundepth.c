//Xander Nuttle
//sundepth.c
//Call: sundepth suns_file sunks_file fastq_file_with_SUNK-containing_reads
//
//This program calculates read depth over SUNs, taking a file detailing SUNs of interest, a file detailing SUNKs overlapping those SUNs, and a fastq file having
//only reads containing SUNKs overlapping the SUNs of interest as inputs. It outputs a file detailing SUNs of interest with a new column added corresponding to
//the read depth over each SUN.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>
#define sunklength 36

void revcomp(char string[]);
int compfun(const void*p1,const void*p2);

//set up SUNK structure
struct sunk
{
	char contig[51];
  long contigloc;
  char seq[sunklength+1];
};

unsigned long bsi(char*kmer,struct sunk*pointer_array,unsigned long num_kmers);

int main(int argc,char*argv[])
{
	//set up SUN structure
	struct sun
  {
    char contig[51];
    long contigloc;
    long alignloc;
		int readdepth;
		int hit;
  };	

	//determine the total number of SUNs
	char dummy[100];
	FILE*sunsfile;
  fpos_t pos;
	long numsuns=0;
	sunsfile=fopen(*(argv+1),"r");
	fgetpos(sunsfile,&pos);
  while((fscanf(sunsfile,"%s %s %s %s %s",dummy,dummy,dummy,dummy,dummy))==5)
  {
    numsuns++;
  }
	fsetpos(sunsfile,&pos);

	//allocate memory to store SUN information
	struct sun*suns;
  long i,j;
	suns=(struct sun*)malloc(numsuns*sizeof(struct sun));
  if(suns==NULL)
  {
    printf("Memory allocation for mip information failed!\n");
    return 1;
  }
	for(i=0;i<numsuns;i++)
	{
		for(j=0;j<51;j++)
		{
			suns[i].contig[j]='\0';
		}
	}

	//read in suns file and store info about each SUN
	i=0;
  while(fscanf(sunsfile,"%s %ld %ld %s %s",&(suns[i].contig),&(suns[i].contigloc),&(suns[i].alignloc),dummy,dummy)==5)
	{
		suns[i].readdepth=0;
		suns[i].hit=0;
		i++;
	}
	fclose(sunsfile);

	//determine total number of SUNKs
	FILE*sunksfile;	
	long numsunks=0;
	sunksfile=fopen(*(argv+2),"r");
	fgetpos(sunksfile,&pos);
	while((fscanf(sunksfile,"%s %s %s",dummy,dummy,dummy))==3)
  {
    numsunks++;
  }
	fsetpos(sunksfile,&pos);

	//allocate memory to store SUNK information
	struct sunk*sunks;
	sunks=(struct sunk*)malloc((2*numsunks+2)*sizeof(struct sunk));
  if(sunks==NULL)
  {
    printf("Memory allocation for mip information failed!\n");
    return 1;
  }
	for(i=0;i<(2*numsunks+2);i++)
  {
    for(j=0;j<51;j++)
    {
      sunks[i].contig[j]='\0';
    }
  }

	//read in sunks file and store info about each SUNK
	char firststring[]="!";
  char laststring[]="~";
	strncpy(sunks[0].contig,firststring,sunklength);
	strncpy(sunks[0].seq,firststring,sunklength);
	sunks[0].contigloc=0;
	i=1;
  while(fscanf(sunksfile,"%s %ld %s",&(sunks[i].contig),&(sunks[i].contigloc),&(sunks[i].seq))==3)
  {
    strncpy(sunks[i+1].contig,sunks[i].contig,50);
		sunks[i+1].contigloc=sunks[i].contigloc;
		strncpy(sunks[i+1].seq,sunks[i].seq,sunklength);
		revcomp(sunks[i+1].seq);
		i+=2;
  }
  fclose(sunksfile);
	strncpy(sunks[i].contig,laststring,sunklength);
  strncpy(sunks[i].seq,laststring,sunklength);
  sunks[i].contigloc=0;
	qsort(sunks,(2*numsunks+2),sizeof(struct sunk),compfun);

	//process fastq file and calculate read depth over all SUNs
	FILE*fastq;
	int readlength;
  unsigned long loc;
  char line1[2001];
  char line2[2001];
  char line3[2001];
  char line4[2001];
  fastq=fopen(*(argv+3),"r");
	fgetpos(fastq,&pos);
  fgets(line2,2000,fastq);
  for(i=0;i<2001;i++)
  {
  	line2[i]='\0';
 	}
  fgets(line2,2000,fastq);
  fsetpos(fastq,&pos);
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
  while(fgets(line1,2000,fastq))
  {
  	for(j=0;j<numsuns;j++)
		{
			suns[j].hit=0;
		}
		fgets(line2,2000,fastq);
  	fgets(line3,2000,fastq);
  	fgets(line4,2000,fastq);
   	for(i=0;i<(readlength-sunklength+1);i++)
    {
    	strncpy(kmer,(line2+i),sunklength);
      loc=bsi(kmer,sunks,(2*numsunks+2));
			if(loc!=ULONG_MAX)
      {
				for(j=0;j<numsuns;j++)
    		{
					if((strcmp(sunks[loc].contig,suns[j].contig)==0)&&(sunks[loc].contigloc>=(suns[j].contigloc-(sunklength-1)))&&(sunks[loc].contigloc<=(suns[j].contigloc)))
      		{
        		suns[j].hit=1;
      		}
    		}      	
      }
    }
		for(j=0;j<numsuns;j++)
		{
			suns[j].readdepth+=suns[j].hit;
		}
    if(feof(fastq))
    {
    	break;
    }
  }
  fclose(fastq);

	//print to output data for all SUNs
	FILE*out;
  out=fopen("brkpt.suns.depth","w");
	for(j=0;j<numsuns;j++)
	{
		fprintf(out,"%s\t%ld\t%ld\t%d\n",suns[j].contig,suns[j].contigloc,suns[j].alignloc,suns[j].readdepth);
	}
  
	//clean up and exit
	fclose(out);
	free(suns);
	free(sunks);
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
	const struct sunk*a1=p1;
  const struct sunk*a2=p2;
  return(strncmp(a1->seq,a2->seq,sunklength));
}

unsigned long bsi(char*kmer,struct sunk*pointer_array,unsigned long num_kmers)
{
  unsigned long start=0;
  unsigned long end=num_kmers;
  unsigned long index;
  int i;
  while(start<end)
  {
    index=(start+end)/2;
    i=strncmp(kmer,pointer_array[index].seq,sunklength);
    if(i==0) {return index;}
    else if(i>0) {start=index+1;}
    else {end=index-1;}
  }
  index=(start+end)/2;
  i=strncmp(kmer,pointer_array[index].seq,sunklength);
  if(i==0) {return index;}
  return ULONG_MAX;
}
