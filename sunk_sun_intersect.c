//Xander Nuttle
//sunk_sun_intersect.c
//Call: sunk_sun_intersect suns_file sunks_file
//
//Taking as inputs a file detailing SUNs of interest and another file detailing SUNKs over regions of interest, this program outputs to a new file all SUNKs
//overlapping a SUN of interest.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>

void revcomp(char string[]);

int main(int argc,char*argv[])
{
	//set up SUN structure
	struct sun
  {
    char contig[51];
    long contigloc;
    long alignloc;
		int readdepth;
  };	

	//determine the total number of SUNs
	char dummy[50];
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

	//allocate memory to store sun information
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
		i++;
	}
	fclose(sunsfile);

	//determine SUNK size and total number of SUNKs
	FILE*sunksfile;	
	long numsunks=0;
	char kmer[201];
	int sunklength;
	for(i=0;i<201;i++)
	{
		kmer[i]='\0';
	}
	sunksfile=fopen(*(argv+2),"r");
	fgetpos(sunksfile,&pos);
	while((fscanf(sunksfile,"%s %s %s",dummy,dummy,kmer))==3)
  {
    numsunks++;
  }
	fsetpos(sunksfile,&pos);
	sunklength=strlen(kmer);

	//set up SUNK structure
	struct sunk
  {
    char contig[51];
    long contigloc;
    char seq[sunklength+1];
  };

	//allocate memory to store sunk information
	struct sunk*sunks;
	sunks=(struct sunk*)malloc(numsunks*sizeof(struct sunk));
  if(sunks==NULL)
  {
    printf("Memory allocation for mip information failed!\n");
    return 1;
  }
	for(i=0;i<numsunks;i++)
  {
    for(j=0;j<51;j++)
    {
      sunks[i].contig[j]='\0';
    }
  }

	//read in sunks file and store info about each SUNK
	i=0;
  while(fscanf(sunksfile,"%s %ld %s",&(sunks[i].contig),&(sunks[i].contigloc),&(sunks[i].seq))==3)
  {
		i++;
  }
  fclose(sunksfile);

	//set up output file 
	FILE*out;
	char output_base_name[101];
  for(i=0;i<100;i++)
  {
    output_base_name[i]='\0';
  }
  char output_file_extension[10]=".sunsunks";
  strncpy(output_base_name,*(argv+2),strrchr(*(argv+2),'.')-(*(argv+2)));
  strcat(output_base_name,output_file_extension);
  out=fopen(output_base_name,"w");	

	//determine which SUNKs overlap at least 1 SUN and print these SUNKs to an output file
	for(i=0;i<numsunks;i++)
  {
    for(j=0;j<numsuns;j++)
    {
			if((strcmp(sunks[i].contig,suns[j].contig)==0)&&(sunks[i].contigloc>=(suns[j].contigloc-(sunklength-1)))&&(sunks[i].contigloc<=(suns[j].contigloc)))
			{
				fprintf(out,"%s\t%ld\t%s\n",sunks[i].contig,sunks[i].contigloc,sunks[i].seq);
				break;
			}
    }
  }

	//clean up and exit
	fclose(out);
	return 0;
}
