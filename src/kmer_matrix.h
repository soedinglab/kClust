#ifndef CM_KMER_MATRIX_H
#define CM_KMER_MATRIX_H

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h> 

#include "simple_list.h"
#include "my_exception.h"


class Kmer_matrix{
	
	public:
	
	struct entry_score{
		public:
		entry_score(int kmer_idx=0, float score=0.0f): kmer_idx(kmer_idx), score(score){}
		int kmer_idx;
		float score;
	};
	
	Kmer_matrix(int k, float scale);
	
	~Kmer_matrix();
	
	Simple_list<entry_score>* get_kmer_list(int kmer_idx);
	
	void fill_kmer_matrix(std::string file);
	
	int get_aa_dim(){ return AMINOACID_DIM;}
	
	size_t* get_aa_powers() { return aa_powers;}
	
	private:
	
	int k;
	float scale;
	Simple_list<entry_score> ** scores;
	int AMINOACID_DIM;
	int size;
	size_t *aa_powers;
	char *int2aa;
	int *aa2int;

};
#endif