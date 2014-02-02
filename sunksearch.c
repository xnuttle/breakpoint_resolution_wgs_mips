/*
 * sunksearch.c
 *
 *  Created on: Nov 3, 2011
 *      Author: Xander
 */

#include<stdio.h>
#include<ctype.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>

int compfun(const void*p1,const void*p2);
unsigned long bsi(char*kmer,char**pointer_array,unsigned long num_kmers);
void revcomp(char*string);
int mersize=36;

int main(int argc,char*argv[])
{
	
	if(argc<5)
	{
		printf("Not enough input!\n"); 
		return 0;
	}
	FILE*files[argc-1];
	int i;
	for(i=0;i<(argc-2);i++)
	{
	files[i]=fopen((*(argv+1+i)),"r");
	}
	
	int masking=1;
	unsigned long genomesize=3400000000,contigs_size=10000000,maskstart,maskend; //genomesize=3400000000 is good
	char*genome,*contiggenome,*contiggenomecopy;
	char**pointers;
	char contignames[argc-4][51];
	genome=(char*)malloc(genomesize*sizeof(char));
	contiggenome=(char*)malloc(contigs_size*sizeof(char));
	contiggenomecopy=(char*)malloc(contigs_size*sizeof(char));
	pointers=(char**)malloc((genomesize+contigs_size)*sizeof(char*));
	if((genome==NULL)||(pointers==NULL)||(contiggenome==NULL)||(contiggenomecopy==NULL))
	{
		printf("Memory allocation failed!\n");
		exit(EXIT_FAILURE);
	}
	int r=-1;
	char ch;
	char chr[31],maskchr[31],string[mersize+1],stringRC[mersize+1],contig[51];
	unsigned long k=1,m=0,n=1,p=0,q=0,s=0,loc,locRC;
	fscanf(files[argc-3],"%s%lu%lu",maskchr,&maskstart,&maskend);
	char firststring[]="!";
	pointers[0]=firststring;
	while((ch=getc(files[0]))!=EOF)
	{
		k=1;
		if(ch=='>')
		{
			fscanf(files[0],"%s",&chr);
			genome[m]='X';
			pointers[m+1]=(genome+m);	
			m++;
			while((strcmp(chr,maskchr))==0)
			{
				while(k<maskstart)
				{
					ch=getc(files[0]);
					if(isalpha(ch))
					{
						genome[m]=ch;
						pointers[m+1]=(genome+m);
						m++;
						k++;
					}
				}
				genome[m]='M';
				pointers[m+1]=(genome+m);
				m++;
				while(k<=maskend)
				{
					ch=getc(files[0]);
					if(isalpha(ch))
					{
						k++;
					}
				}
				masking=fscanf(files[argc-3],"%s %lu %lu",maskchr,&maskstart,&maskend);
				if(masking<3)
				{
					ch=getc(files[0]);
					break;
				}
			}
		}	
		if(isalpha(ch))
		{
			genome[m]=ch;
			pointers[m+1]=(genome+m);
			m++;
		}
	}
	genome[m]='\0';
	
	while(n<(argc-3))
	{
		while((ch=getc(files[n]))!=EOF)
		{
			if(ch=='>')
			{
				fscanf(files[n],"%s",contignames[n-1]);
				contiggenome[p]='X';
				pointers[m+1]=(contiggenome+p);	
				contiggenomecopy[p]=contiggenome[p];
				p++;
				m++;
				ch=getc(files[n]);
			}
			if(isalpha(ch))
			{
				contiggenome[p]=ch;
				contiggenomecopy[p]=contiggenome[p];
				pointers[m+1]=(contiggenome+p);
				p++;
				m++;
			}
		}
		n++;
	}

	qsort(pointers,m,sizeof(char*),compfun);

	files[i]=fopen((*(argv+1+i)),"w");	
	
	for(q=0;q<p-(mersize-1);q++)
	{
		if(*(contiggenomecopy+q)=='X')
		{
			r++;
			s=0;
		}
		strncpy(string,contiggenomecopy+q,mersize);
		strcpy(stringRC,string);
		revcomp(stringRC);
		if(strpbrk(string,"XN")==NULL)
		{
			loc=bsi(string,pointers,m);
			locRC=bsi(stringRC,pointers,m);
			if((locRC==ULONG_MAX)&&(strncmp(string,pointers[loc-1],mersize)!=0)&&(strncmp(string,pointers[loc+1],mersize)!=0))
			{
				fprintf(files[argc-2],"%s\t%lu\t%s\n",contignames[r],s,string);
			}
		}
		s++;
	}

	for(i=0;i<(argc-1);i++)
	{
	fclose(files[i]);
	}
	free(genome);
	free(contiggenome);
	free(contiggenomecopy);
	free(pointers);
	return 0;
}









int compfun(const void*p1,const void*p2)
{
	char*const*a1=p1;
	char*const*a2=p2;
	char s1[mersize+1],s2[mersize+1];
	strncpy(s1,*a1,mersize);
	strncpy(s2,*a2,mersize);
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
		i=strncmp(kmer,pointer_array[index],mersize);
		if(i==0) {return index;}
		else if(i>0) {start=index+1;}
		else {end=index-1;}
	}
	index=(start+end)/2;
	i=strncmp(kmer,pointer_array[index],mersize);
	if(i==0) {return index;}
	return ULONG_MAX;
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
		case 'A': 	base='T';
					break;
		case 'C': 	base='G';
					break;
		case 'G':	base='C';
					break;
		case 'T':	base='A';
					break;
		}
		newstring[c]=base;
	}
	newstring[a]='\0';
	strcpy(string,newstring);
}
