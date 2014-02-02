//Xander Nuttle
//calculate_target_GC.c
//Call: calculate_target_GC mip.sequences
//
//Calculates and outputs the GC content of putative target regions, taking an input file having the sequences of each of these regions each on a separate line.

#include<stdio.h>

int main(int argc,char*argv[])
{
	FILE*in;
	double numbases;
	int numGC;
	char base;
	in=fopen(*(argv+1),"r");
	while((base=getc(in))!=EOF)
	{
		numGC=0;
		numbases=0.0;
		while(base!='\n')
		{
			numbases++;
			if((base=='A')||(base=='T'))
				;
			else
				numGC++;
			base=getc(in);
		}
		printf("%lf\n",(double)(numGC/numbases));
	}
	fclose(in);
	return 0;
}
