//Xander Nuttle
//pickmips.c
//Given a set of potential MIP targets to pick from and a set of scored SUNs for a given paralog, this program selects MIPs to order.
//Call: ./pickmips mipdesign_file SUNs_scored_file (int)number_of_paralogs_in_gene_family
//
//Compiling command: gcc pickmips.c -o pickmips -lm $CPPFLAGS $LDFLAGS
//
//The program systematically and aggressively selects MIPs to order. It attempts to pick the maximum number of MIPs that never target the same SUN
//to cover as many SUNs as possible while prioritizing MIPs having high paralog-specificity. Several SUNs get a score of 0 simply because they are
//located within masked sequence. However, most such SUNs are still useful for paralog-specific copy number genotyping. The program implements a
//greedy algorithm to attempt to find the maximum set packing, treating the universe U as the set of all SUNs and each MIP as a set S containing
//targeted SUNs as elements. The goal is then to select the maximum number of MIPs (sets S) without allowing any selected MIP (set) to target the same
//SUN (to share a common element). The reason for this is that separate MIP data points can be treated as independent as long as no targeted SUNs are
//shared; the only autocorrelation in the data will be due to a paralog-specific copy number state remaining constant across a significant genomic
//spatial extent. By considering at first only MIPs that allow all paralogs to be distinguished, attempting to find the maximum set packing, and repeating
//this process until finsishing with MIPs that allow only 1 paralog to be distinguished from all others, MIPs with high paralog-specificity are effectively
//prioritized. Finally, there are some theoretical cases where this approach could result in many SUNs not being covered due to the requirement of no
//shared targeted SUNs in the final MIP set. Thus, the program reports the number of coverable SUNs (if all MIPs are selected), the number of SUNs
//actually covered by the selected MIP set, and information on coverable SUNs that were not covered. This allows for manual refinement of the selected
//MIP set if too many SUNs or attractive high-scoring SUNs end up not targeted by the MIP set selected by the program.
//
//Note: the MIP design file should have 21 columns, including 1 with GC percentage, 1 with the MIP oligo sequence with backbone, and none with SNP info.
//Note: the program works for a maximum MIP arm length of 24 bases; it could be generalized or changed if need be to account for longer arms

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>
#include<math.h>

int main(int argc,char*argv[])
{
	//determine MIP length, MIP target length, and number of MIPs to be considered for selection
	FILE*designfile;
	fpos_t start_designfile;
	int miplength;
	int targetlength;
	char dummy[301];
	long i;
	long nummips=0;
	designfile=fopen(*(argv+1),"r");
	fgetpos(designfile,&start_designfile);
	fscanf(designfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s",dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy);
	for(i=0;i<301;i++)
  {
    dummy[i]='\0';
  }
	fscanf(designfile,"%s",dummy);
	targetlength=strlen(dummy);
	fscanf(designfile,"%s %s %s %s %s %s",dummy,dummy,dummy,dummy,dummy,dummy);
	for(i=0;i<301;i++)
  {
    dummy[i]='\0';
  }
  fscanf(designfile,"%s",dummy);
	miplength=strlen(dummy);
	fsetpos(designfile,&start_designfile);
	while(fscanf(designfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)==21)
	{
		nummips++;
	}
	fsetpos(designfile,&start_designfile);

	//set up MIP structure
	int max_armlength=24;
	struct mip
	{
		long num;
		int score;
		char master_target_name[51];
		long ext_arm_start;
		long ext_arm_end;
		char ext_arm_seq[max_armlength+1];
		int ext_arm_cn;
		long lig_arm_start;
		long lig_arm_end;
		char lig_arm_seq[max_armlength+1];
		int lig_arm_cn;
		long target_start;
		long target_end;
		char target_seq[targetlength+1];
		long feature_start;
		long feature_end;
		double feature_count;
		char strand;
		double target_gc;
		char note[51];
		char mipseq[miplength+1];
		int mip_picked;
		long num_plogs_dist;
		int num_ovlp_samespec;
		int num_ovlp_lowspec;
		double sunscore;
	};

	//allocate memory to store MIP information
	struct mip*mips;
	mips=(struct mip*)malloc(nummips*sizeof(struct mip));
	if(mips==NULL)
	{
		printf("Memory allocation for mip information failed!\n");
		return 1;
	}

	//read in MIP information from MIP design file
	i=0;
	while(fscanf(designfile,"%ld %d %s %ld %ld %s %d %ld %ld %s %d %ld %ld %s %ld %ld %lf %c %lf %s %s",&(mips[i].num),&(mips[i].score),mips[i].master_target_name,&(mips[i].ext_arm_start),&(mips[i].ext_arm_end),mips[i].ext_arm_seq,&(mips[i].ext_arm_cn),&(mips[i].lig_arm_start),&(mips[i].lig_arm_end),mips[i].lig_arm_seq,&(mips[i].lig_arm_cn),&(mips[i].target_start),&(mips[i].target_end),mips[i].target_seq,&(mips[i].feature_start),&(mips[i].feature_end),&(mips[i].feature_count),&(mips[i].strand),&(mips[i].target_gc),mips[i].note,mips[i].mipseq)==21)
	{
		mips[i].mip_picked=0;
		mips[i].num_plogs_dist=0;
		mips[i].num_ovlp_samespec=0;
		mips[i].num_ovlp_lowspec=0;
		mips[i].sunscore=0.0;
		i++;
	}
	fclose(designfile);
	
	//determine number of SUNs
	FILE*sunfile;
  fpos_t start_sunfile;
  long numsuns=0;
  sunfile=fopen(*(argv+2),"r");
  fgetpos(sunfile,&start_sunfile);
  while(fscanf(sunfile,"%s %s %s %s %s",dummy,dummy,dummy,dummy,dummy)==5)
  {
    numsuns++;
  }
  fsetpos(sunfile,&start_sunfile);	

	//set up SUN structure
	struct sun
	{
		char contig_name[51];
		long contigloc;
		long alignloc;
		double score;
		char note[51];
		int sun_coverable;
		int sun_covered;
	};

	//allocate memory to store SUN information
	struct sun*suns;
	suns=(struct sun*)malloc(numsuns*sizeof(struct sun));
	if(suns==NULL)
  {
    printf("Memory allocation for sun information failed!\n");
    return 1;
  }

	//read in SUN information from scored SUNs file
	i=0;
	while(fscanf(sunfile,"%s %ld %ld %lf %s",suns[i].contig_name,&(suns[i].contigloc),&(suns[i].alignloc),&(suns[i].score),suns[i].note)==5)
	{
		suns[i].sun_coverable=0;
		suns[i].sun_covered=0;
		i++;
	}
	fclose(sunfile);

	//calculate maximum SUN score
	double max_sunscore=0.0;
	long j;
	for(j=0;j<numsuns;j++)
	{
		if(suns[j].score>max_sunscore)
		{
			max_sunscore=suns[j].score;
		}
	}

	//get information from the command line regarding the number of paralogs in the gene family
	long num_paralogs;
	num_paralogs=strtol(*(argv+3),NULL,10);

	//create and initialize coverage matrix to show MIP-SUN targeting; MIPs are rows, SUNs are columns
	int(*coverage)[numsuns]=malloc(nummips*sizeof(*coverage)); 
	if(coverage==NULL)
	{
		printf("Memory allocation for sun information failed!\n");
		return 1;
	}
	for(i=0;i<nummips;i++)
	{
		for(j=0;j<numsuns;j++)
		{
			coverage[i][j]=0;
		}
	}

	//go through all MIPs to fill out coverage matrix: entry becomes 1 if MIP i targets SUN j
	//also determine paralog-specificity of each MIP i and coverability of each SUN j
	char diff_contig_names[targetlength*num_paralogs][51];
	long diff_contig_locs[targetlength*num_paralogs];
	long k,m,num_contigs_seen;
	int contig_seen;
	for(i=0;i<nummips;i++) //i is MIP
	{
		num_contigs_seen=0;
		for(k=0;k<targetlength;k++)
		{
			for(m=0;m<51;m++)
			{
				diff_contig_names[k][m]='\0';
			}
			diff_contig_locs[k]=LONG_MIN;
		}
		k=0;
		for(j=0;j<numsuns;j++) //j is SUN
		{
			if((suns[j].alignloc>=mips[i].target_start)&&(suns[j].alignloc<=mips[i].target_end))
			{
				coverage[i][j]=1;
				suns[j].sun_coverable=1;
				if(strncmp(suns[j].note,"SUN_is_masked",13)==0) //many of these SUNs are likely to be real and near fixed for a particular paralog, so their scores should not be zero
				{
					mips[i].sunscore+=((2.0/3.0)*max_sunscore); //arbitrarily give masked SUNs scores of (2/3) times the maximum observed SUN score
				}
				else if(strncmp(suns[j].note,"no_SUNK_ovlp",12)==0) //many of these SUNs are likely part of longer SUNKs even though they are not part of 30-mer SUNKs, so their scores should not be zero
				{
					mips[i].sunscore+=(0.5*max_sunscore); //arbitrarily give SUNs not part of 30-mer SUNKs SUN scores of (1/2) times the maximum observed SUN score
				}
				else
				{
					mips[i].sunscore+=suns[j].score;
				}
				if(k==0) //k is SUN targeted by MIP
				{
					strncpy(diff_contig_names[k],suns[j].contig_name,50);
					diff_contig_locs[k]=suns[j].contigloc;
					num_contigs_seen++;
				}
				else
				{
					contig_seen=0;
					for(m=0;m<num_contigs_seen;m++) //m is distinct contig
					{
						if((strncmp(suns[j].contig_name,diff_contig_names[m],50)==0)&&((suns[j].contigloc-diff_contig_locs[m])<=(targetlength-1))&&((suns[j].contigloc-diff_contig_locs[m])>=(1-targetlength)))
						{
							contig_seen=1;
						}
					}
					if(!(contig_seen))
					{
						strncpy(diff_contig_names[num_contigs_seen],suns[j].contig_name,50);
						diff_contig_locs[num_contigs_seen]=suns[j].contigloc;
						num_contigs_seen++;
					}
				}
				k++;
			}
		}
		if(num_contigs_seen>=(num_paralogs-1))
		{
			mips[i].num_plogs_dist=num_paralogs;
		}
		else
		{
			mips[i].num_plogs_dist=num_contigs_seen;
		}
	}

	//go through set of all MIPs and for each distinct SUN set select a best MIP and eliminate all other MIPs targeting that SUN set
	long mips_remaining,bestmip;
	int mip_is_equiv;
	for(i=0;i<nummips;i++)
	{
		if(mips[i].mip_picked==0)
		{
			bestmip=i;
			mips_remaining++;
			for(k=i+1;k<nummips;k++)
			{
				mip_is_equiv=1;
				for(j=0;j<numsuns;j++)
				{
					if(coverage[i][j]!=coverage[k][j])
					{
						mip_is_equiv=0;
						break; //MIPs i and k target distinct SUN sets
					}
				}
				if(!mip_is_equiv)
				{
					break;
				}
				else
				{
					if(mips[k].score>mips[bestmip].score) //first criterion for picking the best MIP from a set of MIPs targeting the same SUN set: MIP design score
					{
						mips[bestmip].mip_picked=-1;
						bestmip=k;
					}
					else if(mips[k].score<mips[bestmip].score)
					{
						mips[k].mip_picked=-1;
					}
					else
					{
						if((fabs(mips[k].target_gc-0.45))<(fabs(mips[bestmip].target_gc-0.45))) //second criterion: MIP target GC content
						{
            	mips[bestmip].mip_picked=-1;
							bestmip=k;
						}
						else if((fabs(mips[k].target_gc-0.45))>(fabs(mips[bestmip].target_gc-0.45)))
						{
							mips[k].mip_picked=-1;
						}
						else
						{
							if((mips[k].ext_arm_cn+mips[k].lig_arm_cn)<(mips[bestmip].ext_arm_cn+mips[bestmip].lig_arm_cn)) //third criterion: MIP arm total copy count
							{
								mips[bestmip].mip_picked=-1;
								bestmip=k;
							}
							else
							{
								mips[k].mip_picked=-1;
							}
						}
					}
				}	
			}
			i=k-1;
		}
	}

	//pick MIPs
	long max_specificity;
	int mip_ovlp;
	while(mips_remaining)
	{
		//determine maximum paralog-specificity of remaining MIPs
		max_specificity=0;
		for(i=0;i<nummips;i++)
		{
			if(((mips[i].mip_picked==0)||(mips[i].mip_picked==2))&&(mips[i].num_plogs_dist>max_specificity))
			{
				max_specificity=mips[i].num_plogs_dist;
			}
		}
		//for each remaining MIP having this maximum paralog-specificity, calculate the number of remaining MIPs with the same specificity it intersects
		//also calculate the number of remaining MIPs with lower paralog-specificity it intersects
		//if it does not intersect any remaining MIPs with the same specificity, pick it and eliminate any remaining MIPs it intersects; otherwise flag it
		for(i=0;i<nummips;i++)
		{
			if((mips[i].mip_picked==0)&&(mips[i].num_plogs_dist==max_specificity))
			{
				for(k=0;k<nummips;k++)
				{
					if((i!=k)&&((mips[k].mip_picked==0)||(mips[k].mip_picked==2)))
					{
						for(j=0;j<numsuns;j++)
						{
							if((coverage[i][j]==coverage[k][j])&&(coverage[i][j]))
							{
								if(mips[k].num_plogs_dist==max_specificity)
								{
									mips[i].num_ovlp_samespec++;
									break;
								}
								else
								{
									mips[i].num_ovlp_lowspec++;
									break;
								}
							}
						}
					}
				}
				if(mips[i].num_ovlp_samespec==0)
				{
					mips[i].mip_picked=1;
					mips_remaining--;
					for(j=0;j<numsuns;j++)
					{
						if(coverage[i][j])
						{
							suns[j].sun_covered=1;
						}
					}
					for(k=0;k<nummips;k++)
					{
						if(mips[k].mip_picked==0)
						{
							for(j=0;j<numsuns;j++)
            	{
              	if((coverage[i][j]==coverage[k][j])&&(coverage[i][j]))
              	{
									mips[k].mip_picked=-1;
									mips_remaining--;
									break;
								}
							}
						}
					}
				}
				else
				{
					mips[i].mip_picked=2;
				}
			}
		}
		//resolve conflicts between intersecting MIPs having the same specificity
		for(i=0;i<nummips;i++)
		{
			if(mips[i].mip_picked==2)
			{
				bestmip=i;
				for(k=i+1;k<nummips;k++)
      	{
        	if(mips[k].mip_picked==2)
					{
						mip_ovlp=0;
						for(j=0;j<numsuns;j++)
        		{
          		if((coverage[bestmip][j]==coverage[k][j])&&(coverage[bestmip][j]))
          		{
            		mip_ovlp=1;
            		break; //MIPs bestmip and k target at least 1 common SUN
          		}
        		}
        		if(!mip_ovlp)
        		{
          		break;
        		}
        		else
        		{
							if(mips[k].num_ovlp_samespec<mips[bestmip].num_ovlp_samespec) //first criterion for picking the best MIP from a set of MIPs targeting at least 1 common SUN: number of remaining conflicting MIPs with the same specificity
							{
								mips[bestmip].mip_picked=-1;
								bestmip=k;
								mips_remaining--;			
							}
							else if(mips[k].num_ovlp_samespec>mips[bestmip].num_ovlp_samespec)
							{
								mips[k].mip_picked=-1;
								mips_remaining--;
							}
							else
							{
								if(mips[k].num_ovlp_lowspec<mips[bestmip].num_ovlp_lowspec) //second criterion for picking the best MIP from a set of MIPs targeting at least 1 common SUN: number of remaining conflicting MIPs
								{
									mips[bestmip].mip_picked=-1;
									bestmip=k;
									mips_remaining--;
								}
								else if(mips[k].num_ovlp_lowspec>mips[bestmip].num_ovlp_lowspec)
								{
									mips[k].mip_picked=-1;
									mips_remaining--;
								}
								else
								{
									if(mips[k].sunscore>mips[bestmip].sunscore) //third criterion for picking the best MIP from a set of MIPs targeting at least 1 common SUN: combined SUN score
									{
										mips[bestmip].mip_picked=-1;
										bestmip=k;
										mips_remaining--;
									}
									else
									{
										mips[k].mip_picked=-1;
                  	mips_remaining--;
									}
								}
							}
						}
					}
				}
				mips[bestmip].mip_picked=1;
				mips_remaining--;
				for(j=0;j<numsuns;j++)
        {
        	if(coverage[bestmip][j])
          {
          	suns[j].sun_covered=1;
          }
        }
        for(k=0;k<nummips;k++)
        {
        	if(mips[k].mip_picked==0)
          {
          	for(j=0;j<numsuns;j++)
            {
            	if((coverage[bestmip][j]==coverage[k][j])&&(coverage[bestmip][j]))
              {
              	mips[k].mip_picked=-1;
                mips_remaining--;
								break;
              }
            }
          }
        }			
				i=k-1;
			}	
		}
	}

	//print selected MIPs to output file
	FILE*picksfile;
	char output_base_name[101];
	for(i=0;i<100;i++)
	{
		output_base_name[i]='\0';
	}
  char output_file_extension[10]=".mippicks";
	strncpy(output_base_name,*(argv+1),strrchr(*(argv+1),'.')-(*(argv+1)));
  strcat(output_base_name,output_file_extension);
  picksfile=fopen(output_base_name,"w");
	for(i=0;i<nummips;i++)
	{
		if(mips[i].mip_picked==1)
		{
			fprintf(picksfile,"%ld\t%d\t%s\t%ld\t%ld\t%s\t%d\t%ld\t%ld\t%s\t%d\t%ld\t%ld\t%s\t%ld\t%ld\t%.1lf\t%c\t%lf\t%s\t%s\n",mips[i].num,mips[i].score,mips[i].master_target_name,mips[i].ext_arm_start,mips[i].ext_arm_end,mips[i].ext_arm_seq,mips[i].ext_arm_cn,mips[i].lig_arm_start,mips[i].lig_arm_end,mips[i].lig_arm_seq,mips[i].lig_arm_cn,mips[i].target_start,mips[i].target_end,mips[i].target_seq,mips[i].feature_start,mips[i].feature_end,mips[i].feature_count,mips[i].strand,mips[i].target_gc,mips[i].note,mips[i].mipseq);
		}
	}
	fclose(picksfile);
	
	//print information on SUN coverage to screen (standard output)
	long num_suns_coverable=0,num_suns_covered=0;
	for(j=0;j<numsuns;j++)
	{
		if(suns[j].sun_coverable)
		{
			num_suns_coverable++;
		}
		if(suns[j].sun_covered)
		{
			num_suns_covered++;
		}
	}
	printf("Total_number_of_SUNs: %ld\nNumber_of_SUNs_coverable: %ld\nNumber_of_SUNs_covered: %ld (%lf%)\n\nCoverable_SUNs_not_covered:\n",numsuns,num_suns_coverable,num_suns_covered,((((double)num_suns_covered/(double)num_suns_coverable)*100.0)));
	for(j=0;j<numsuns;j++)
	{
		if((suns[j].sun_coverable)&&(!(suns[j].sun_covered)))
		{
			printf("%s\t%ld\t%ld\t%lf\t%s\n",suns[j].contig_name,suns[j].contigloc,suns[j].alignloc,suns[j].score,suns[j].note);
		}
	}

	//clean up and exit program
	free(mips);
	free(suns);
	free(coverage);
	return 0;
}
