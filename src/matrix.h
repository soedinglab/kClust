#ifndef CM_MATRIX_H
#define CM_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Amino acid substitution matrix class
//Exact values for substitution scores are computed from the raw BLOSUM data
//Export BLOSUM_MATRICES='/dir/' to point to the directory where the blosumXX.out files reside
//The files are available at ftp://ftp.ncbi.nih.gov/repository/blocks/unix/blosum/blosum.tar.Z
//The dummy 'static_blosum62' can be alternatively used

#include <sstream>      //stringstream
#include <cmath>        //log, log10
#include <fstream>      //file streams
#include <iostream>     //ostream, istream
#include <cstdio>       //printf, fprintf
#include <cstdlib>      //getenv

#include "sequence.h"
#include "wnconj.h"
#include "wnmem.h"
#include "dummy_matrix.h"

#define CONJ_GRAD_ITERATIONS 100
#define CONJ_GRAD_NORM 1e-10

class Sequence;

class Matrix{
	public:
		struct parameter_wrapper{
			float** matrix;
			float* p_query;
		};
		
		struct dimer{
			unsigned int idx;
			float score;
		};
		

		enum mtype{
			blosum30  = 0,
			blosum35  = 1,
			blosum40  = 2,
			blosum45  = 3,
			blosum50  = 4,
			blosum55  = 5,
			blosum60  = 6,
			blosum62  = 7,
			blosum65  = 8,
			blosum70  = 9,
			blosum75  = 10,
			blosum80  = 11,
			blosum85  = 12,
			blosum90  = 13,
			blosum95  = 14,
			blosum100 = 15,
			static_blosum62   = 16
		};
		Matrix(const mtype t=Matrix::static_blosum62) throw (std::exception);
		~Matrix();	

		//calls modified conjugate gradient method of wnlib (void* added)
		void wn_rescale(Sequence*) throw (std::exception);

		//rounds values in matrix to integer values
		void round_bit_scores();

		//pointers to data 
		const float** const get_matrix() const;
		const float* const get_diagonal_scores() const;
		const float** const get_prob_matrix() const;
		const float* const  get_p_query() const;
		const float* const  get_p_back() const;
		const char * const get_int2aa() const;
		const int  * const get_aa2int() const;
		const size_t get_aa_dim() const;

		//initializes p_query (background frequencies for a certain sequence) for sequence s 
		void init_p_query(Sequence *s);
		//prints float values of matrix
		std::ostream& print_p_query(std::ostream&);

		//prints float values of matrix
		std::ostream& print_p_back(std::ostream&);
		//sets matrix to original blosum bit-scores
		void reset();

		//prints float values of matrix
		std::ostream& print(std::ostream&);

		//prints float values of matrix
		std::ostream& print_frequencies(std::ostream&);
		//prints rounded bit-scores of matrix[][]
		std::ostream& print_int(std::ostream&);

		//set matrix[][] with arbitrary values
		void read_matrix(const float src[][21]);
	
		//returns the number of negative factors after rescaling wth conjugate gradient method
		int get_negatives(){return negatives;}

		//returns the squared error of minimization with conjugate gradient method
		double get_error_sum(){return error_sum;}
		
		size_t * get_aa_powers(){return aa_powers;}
		
		Matrix::dimer ** get_dimers() {return _dimers;}
		
		float get_scale() { return scale;}
		
		inline float _log2 (float x) { return log10(x)/0.301029996; }
		
	private:
		const size_t AMINOACID_DIM;
		float **matrix;
		float *diagonal_scores;
		size_t *aa_powers;	
		float **original;
		float *p_background;
		float *p_query;	
		float scale;
		char *int2aa;
		int *aa2int;
		dimer **_dimers;
		int negatives;
		double error_sum;
		
		void _read_blosum_matrix(const char*) throw (std::exception);
		void _copy_dummy();
		void _copy();
		std::ostream& _print_debug(std::ostream&);
		const std::string _get_fn(const mtype); 
		void _init_main_diagonal_scores();
		//methods for the initialisation of presorted dimers
		void _init_dimers();
		dimer** _merge_sort_dimers(dimer **dimers, size_t len);
		
};
//function to minimize for wnlib
double _func( double*, void*); 
//gradient of the function to minimize for wnlib
void _grad( double*, double* , void*);

#endif
