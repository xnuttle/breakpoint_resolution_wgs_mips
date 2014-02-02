#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <string.h>
using namespace std;

//#include<ext/hash_map>
//using namespace __gnu_cxx;
#include <tr1/unordered_map>

//take in a list of FASTA files to analyze
//open up Ian's file of possible MIPs
//add each of the possible MIPs into a hash table
//for each FASTA in the file list, read FASTA, write forward and reverse complement, increment hash
//iterate through ian's file again, calling the hash_table, writing a new file w/ the copy number result

void load_fasta_list(char * file_name, vector<string> &file_list){
	FILE * fin;
	fin = fopen(file_name, "rt");
	char temp_file[512];
	char * temp_file2;
	while (!feof(fin)){
		fgets(temp_file, 512, fin);	
		if (!feof(fin)){
			temp_file2 = strtok(temp_file, "\n");
			file_list.push_back(temp_file2);
		}
	}
	cout<<file_list.size()<<" FASTA files to be analyzed."<<endl;	
}

char reverse_complement(char start, map<char,char> * complements_in){
	if (complements_in == NULL)
	{
		map<char,char> complements;
		complements['A'] = 'T';
		complements['T'] = 'A';
		complements['G'] = 'C';
		complements['C'] = 'G';
		complements['U'] = 'A';
		complements['M'] = 'K';
		complements['R'] = 'Y';
		complements['W'] = 'W';
		complements['S'] = 'S';
		complements['Y'] = 'R';
		complements['K'] = 'M';
		complements['V'] = 'B';
		complements['H'] = 'D';
		complements['D'] = 'H';
		complements['B'] = 'V';
		complements['N'] = 'N';
		complements['-'] = '-';
		return complements[start];
	}
	else{
		return complements_in->at(start);
	}
}

void load_fasta (FILE * fasta, char * genome, char * genomerc, map<char, char> * complements){
	char fasta_header[512];
	char temp;
	fgets(fasta_header, 512, fasta);
	cout<<fasta_header;
	int genome_location = 0;
	while (!feof(fasta)){
		temp = fgetc(fasta);
		if (temp!='\n' && !feof(fasta)){
			if (temp!='-'){
				temp = char(toupper(temp));
			}
			genome[genome_location] = temp;
			genome_location ++;	
		}
	}
	genome[genome_location] = '\0';
	int genomerc_location = 0;
	for (int i = genome_location -1; i>=0; i--){
		genomerc[genomerc_location] = reverse_complement(genome[i], complements);
		genomerc_location++;
	}
}
void populate_mips_table (char * in_file, tr1::unordered_map<string, int> * mip_copy_counts, int * query_lengths){
	FILE * fin;
	fin = fopen(in_file, "rt");
	fpos_t file_position;
	fgetpos(fin, &file_position);
	char file_line[1024];
	char * mip;
	int counter = 0;
	do{
		fgetpos(fin, &file_position);
		fgets(file_line, 1024, fin);
	}while(file_line[0] == '>' && !feof(fin));
	//reset the file position to the first non-header line
	fsetpos(fin, &file_position);
	//read in the first line separately to populate the query_lenghts array
	fgets(file_line, 1024, fin);
	strtok(file_line, "\t");
	for(int i = 0; i<3; i++){
		strtok(NULL, "\t");
	}
	mip = strtok(NULL, "\t");
	query_lengths[0] = strlen(mip);
	mip_copy_counts->insert(make_pair(mip, 0));
	for (int i= 0; i<2; i++){
		strtok(NULL, "\t");
	}
	mip = strtok(NULL, "\t");
	query_lengths[1] = strlen(mip);
	mip_copy_counts->insert(make_pair(mip, 0));
	while (!feof(fin)){
		counter++;
		if (counter%100000 == 0){
			cout<<"Loaded "<<counter<<" lines."<<endl;
		}
		fgets(file_line, 1024, fin);
		strtok(file_line, "\t");
		for(int i = 0; i<3; i++){
			strtok(NULL, "\t");	
		}
		mip = strtok(NULL, "\t");
		//make sure this line is OK even if multiple mips of same sequence
		if (mip!=NULL){
			mip_copy_counts->insert(make_pair(mip, 0));
		}
		for (int i = 0; i<2; i++){
			strtok(NULL, "\t");
		}
		mip = strtok(NULL, "\t");
		if (mip!=NULL){
			mip_copy_counts->insert(make_pair(mip, 0));
		}
	}
	fclose(fin);
}

void update_mips_file(char * in_file, char * out_file, tr1::unordered_map<string, int> &mip_copy_counts){
	vector<string> mips_header;
	FILE * fin;
	fin = fopen(in_file, "rt");
	fpos_t file_position;	
	fgetpos(fin, &file_position);
	char file_line[1024];
	char * mip;
	char * temp_tok1;
	char * temp_tok2;
	int temp_copy_count;
	do{
		fgetpos(fin, &file_position);
		fgets(file_line, 1024, fin);
		if (file_line[0] =='>'){
			mips_header.push_back(string(file_line));
		}
	}while(file_line[0] == '>' && !feof(fin));
	//reset the file position to the first non-header line
	fsetpos(fin, &file_position);
	FILE * fout;
	fout = fopen(out_file, "wt");	
	for (int i=0; i<(mips_header.size()-1); i++){
		fputs(mips_header[i].c_str(), fout);	
	}
	//write new header line
	fputs(">mip_count\tchr\text_probe_start\text_probe_stop\text_probe_sequence\text_probe_copy_count\tlig_probe_start\tlig_probe_stop\tlig_probe_sequence\tlig_probe_copy_count\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tfeature_start_position\tfeature_stop_position\tfeature_mip_count\tprobe_strand\tnotes\n", fout);
        while (!feof(fin)){
		fgets(file_line, 1024, fin);
		if (!feof(fin)){
			temp_tok1 = strtok(file_line, "\n");
			temp_tok2 = strtok(temp_tok1, "\t");
			fputs(temp_tok2, fout);
			fputc('\t', fout);
                	for(int i = 0; i<3; i++){
				temp_tok2 = strtok(NULL, "\t");
				fputs(temp_tok2, fout);
				fputc('\t', fout);
	                }
			mip = strtok(NULL, "\t");
			fputs(mip, fout);
			fputc('\t', fout);
			temp_copy_count =mip_copy_counts[string(mip)];
			fprintf(fout, "%d\t", temp_copy_count);	 	
                	for (int i = 0; i<2; i++){
                       		temp_tok2 = strtok(NULL, "\t");
				fputs(temp_tok2, fout);
				fputc('\t', fout);
                	}
                	mip = strtok(NULL, "\t");
			fputs(mip, fout);
			fputc('\t', fout);
			temp_copy_count = mip_copy_counts[string(mip)];
			fprintf(fout, "%d\t", temp_copy_count);
			temp_tok2 = strtok(NULL, "\t");
			do{
				fputs(temp_tok2, fout);
				fputc('\t', fout);
				temp_tok2 = strtok(NULL, "\t");
			}while(temp_tok2!=NULL);
			fseek(fout, -1, SEEK_CUR);
			fputc('\n', fout);
		}
	}
	fclose(fin);
	fclose(fout);
}

void do_genome_search(char * genome, tr1::unordered_map<string, int> &mip_copy_counts, int * query_lengths){
	string query;
	bool end = false;
	bool end2 = false;
	int i = 0;
	if (query_lengths[0] == query_lengths[1]){
		do{
			if (genome[i+query_lengths[0]-1]!='\0'){
				query = string(&genome[i], query_lengths[0]);
				if (mip_copy_counts.count(query) >0){
					mip_copy_counts[query] +=1;
				}
				i++;
			}
			else{
				end = true;
			}
		}while(!end);
	}
	else{
		do{
			if (!end){
				if (genome[i+query_lengths[0]-1]!='\0'){
					query = string(&genome[i], query_lengths[0]);
					if (mip_copy_counts.count(query) >0){
						mip_copy_counts[query] +=1;
					}
				}
				else{
					end = true;
				}
			}
			if (!end2){
				if (genome[i+query_lengths[1]-1]!='\0'){
					query = string(&genome[i], query_lengths[1]);
					if(mip_copy_counts.count(query) > 0){
						mip_copy_counts[query] +=1;
					}
				}
				else{
					end2 = true;
				}
			}	
			i++;
		}while(!(end && end2));
	}
} 
int main(int argc, char * argv[])
{
	if (argc!= 4)
	{
		cout<<"Usage: genome_compare fasta_list mips_file out_file"<<endl;
		return 0;
	}

	//create the reverse complement table
	map<char,char> complements;
	complements['A'] = 'T';
	complements['T'] = 'A';
	complements['G'] = 'C';
	complements['C'] = 'G';
	complements['U'] = 'A';
	complements['M'] = 'K';
	complements['R'] = 'Y';
	complements['W'] = 'W';
	complements['S'] = 'S';
	complements['Y'] = 'R';
	complements['K'] = 'M';
	complements['V'] = 'B';
	complements['H'] = 'D';
	complements['D'] = 'H';
	complements['B'] = 'V';
	complements['N'] = 'N';
	complements['-'] = '-';
	
	int query_lengths[2];
	//load the MIPs
	tr1::unordered_map<string, int> mip_copy_counts;
	cout<<"Loading mips..."<<endl;	
	vector<string> mips_header;
	populate_mips_table (argv[2], &mip_copy_counts, query_lengths);
	
	//load the list of FASTA files to be searched
	vector<string> file_list;
	load_fasta_list(argv[1], file_list);
	for (int i=0; i<file_list.size(); i++){
		cout<<"Loading genome from FASTA into memory..."<<endl;
		FILE * fasta;
		unsigned long int fasta_length = 0;
		fasta = fopen(file_list.at(i).c_str(), "rt");
		//get the size of the fasta file
		fseek(fasta, 0, SEEK_END);
		fasta_length = ftell(fasta);
		rewind(fasta);
		//allocate memory for the genome based on fasta size
		char * genome = new char[fasta_length];
		char * genomerc = new char[fasta_length];
		load_fasta(fasta, genome, genomerc, &complements);
		fclose(fasta);
		cout<<"FASTA file loaded successfully"<<endl;
		cout<<"Counting mip occurences in forward strand"<<endl;
		do_genome_search(genome, mip_copy_counts, query_lengths);
		cout<<"Counting mip occurences in reverse strand"<<endl;
		do_genome_search(genomerc, mip_copy_counts, query_lengths);
		delete [] genome;
		delete [] genomerc;
	}
	cout<<"Updating MIPS file..."<<endl;
	update_mips_file(argv[2], argv[3], mip_copy_counts);
	return 0;
}
