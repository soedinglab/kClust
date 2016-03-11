#ifndef CM_KMER_H
#define CM_KMER_H
/***************************************************************************
 *   Copyright (C) 2006 by christian mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//module for generatingsimilar k-mers of protein sequences

#include <cstdlib>
#include <string>
#include <list>

#include "matrix.h"
#include "sequence.h"
#include "kmer_seed_list.h"
#include "simple_hash.h"
#include "kmer_matrix.h"

class Kmer{
	
	friend class Clusters;

	public:
		class diag_pos{
			public:
				diag_pos():pos(0),visited(false){}
				diag_pos(size_t pos):pos(pos),visited(false){}
				size_t pos;
				bool visited;
		};		

		struct k_seed_struct{
			char *kmer;	
			k_seed_struct *next;
			char score;
			unsigned short pos;
		};
		
		struct kmer_score{
			float score;
			int   count;	
		};


		Kmer(const size_t k, const double threshold_per_column, bool score_correction, Kmer_matrix* m_4mer = NULL, Kmer_matrix* m_2mer = NULL);
		~Kmer();

		//methods for generation of similar k-mers
		size_t create_kmer_list(Sequence*);
		// only for sequences (dimers not implemented for profiles)
		size_t create_kmer_list_fast(Sequence*) throw (std::exception);

		//PATTERNS: k=6 XX0XXX0X    k=4 XX0XX
		// fast method using dimers	
		// only for sequences (dimers not implemented for profiles)
		size_t create_kmer_list_fast_spaced(Sequence*) throw (std::exception);
		//method to use with an arbitrary substitution matrix
		size_t create_kmer_list_spaced(Sequence*) throw (std::exception);
		
		size_t create_kmer_list_from_matrix(Sequence*) throw (std::exception);

		size_t create_kmer_list_spaced_from_matrix(Sequence*) throw (std::exception);

		//just put all identities in the k-mer list - benchmark of identical k-mer counting
		int create_idents_list(Sequence*);

		//hash based computation of k-mer score - for benchmark purposes
		//TODO now implemented only for sequences
		void get_kmer_score_with(Sequence *, kmer_score&);

		//reports matches of similar k-mers between two sequences
		//TODO now implemented only for sequences
		std::ostream& print_kmer_matches_with(Sequence *, std::ostream &out);

		//computes the k-mer score and count identically to the index-table 
		//and stores them in the referenced structure
		//TODO now implemented only for sequences
//		void get_kmer_score_spaced_with(Sequence *, kmer_score&);

		//computes the k-mer score and count without k-mer score maximization per query position 
		//and stores them in the referenced structure
		//TODO now implemented only for sequences
		void get_full_kmer_score_with(Sequence *, kmer_score&);

		//computes the maximally scoring diagonal
		size_t get_best_diag(size_t *, const size_t);

		//computes the index of a k-mer
		size_t kmer2index(const int*, int begin, int end) const;
		
		//prints the list of generated k-mers
//		void print_kmer_seed_list();

		//number of generated k-mers
		size_t get_list_length();
		size_t get_kmer_count();
	
		//counts the number of identical k-mers - hash based
		size_t count_idents(Sequence *q, Sequence *s, size_t k);

		//prints all matches of identical k-mers
		std::ostream& print_idents(Sequence *q, Sequence *s, size_t k, std::ostream &out);

		//Iterators for traversing the k-mer list
		Kmer_seed_list::Iterator begin();
		Kmer_seed_list::Iterator end();	

		//prints the k-mer with amino acids from its index 
		std::ostream& print_kmer(size_t, std::ostream&);

		//stores the corresponding k-mer of an index in the specified character buffer
		void index2kmer(size_t, char*);
		
		std::string idx2kmer(int idx);

		//returns the approximate memory usage of the instance
		size_t get_memory_usage();

	private:

		//VARIABLES

		//length of kmer
		const size_t k;
		//threshold-score per column
		const float thr;
		//threshold-score for entire k-mer
		const float k_thr;
		//number of amino acids
		size_t AMINOACID_DIM;
		//conversion from numeric aa to character
		const char* int2aa;
		//conversion from character to numeric aa
		const int* aa2int;
		//pointer to substitution matrix
		const float** aa_matrix;
		// pointers to k-mer matrices
		Kmer_matrix* m_4mer;
		Kmer_matrix* m_2mer;
		
		//array of maximally atteinable score of the current kmer from position j..end - k_thr
		//k_thr is substracted to avoid an additional '+' operation in the rec_kmer_count function	
//		float *smax;

		//copy of current kmer
//		char *kmer;

		//current kmer: aa for sequences, profile positions for profiles
//		int * kmer;

		//number of generated k-mers for the current sequence
		size_t kmer_count;		

		//sorted positions of the k-mer in decreasing order of the amino acid identity scores
//		int *sort;

		//copy of diagonal scores of amino acids (speed-up)
		const float * max_scores;

		//21^0, 21^1, 21^2, ...
		size_t *aa_powers;		

		//explicit variables for speed-up
//		int kmer0, kmer1, kmer2, kmer3, kmer4, kmer5, kmer6;
//		float smax0,  smax1, smax2, smax3, smax4, smax5, smax6;
		
		//specialized list object holding the generated k-mers
		Kmer_seed_list seed_list;

		//pointer to sorted 2-mers
		// index i of dimer -> _dimers[i][j] is the j-th best scored dimer with index _dimers[i][j].idx and score _dimers[i][j].score
		Matrix::dimer **_dimers;
		
		float* score_corr;
		
		bool sc;



		//METHODS
		
		void init_to_sequence(Sequence * s);

		//initializes smax for the current k-mer stored in *kmer
		inline void _init_smax(int* kmer, float* smax, float* corr);
		
		//copy identity scores of amino acid into *max_scores
		//required for matrix-rescaling
//		inline void _init_main_diagonal_scores(const float **matrix);

		//sorts the k-mer
		//the amino acid positions of the original k-mer are stored in *sort
		inline void _naiv_kmer_sort(int* kmer, int* sort);
		
		void _sort (float* corr, int* sort);

		//methods for generating similar k-mers with an arbitrary substitution matrix
		//recursive implementation
		void _rec_add_similar_kmers(const size_t, const float, size_t, size_t, int* kmer, int* sort, float* smax);
		//explicit loop implementations
		void _add_similar_kmers(int* kmer, float* smax, int* sort, const size_t i, float* corr);
		void _add_similar_3mers(const size_t, int kmer0, int kmer1, int kmer2, float smax0, float smax1, float smax2, int* sort, float* corr);
		void _add_similar_4mers(const size_t, int kmer0, int kmer1, int kmer2, int kmer3, float smax0, float smax1, float smax2, float smax3, int* sort, float* corr);
		void _add_similar_5mers(const size_t, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, float smax0, float smax1, float smax2, float smax3, float smax4, int* sort, float* corr);
		void _add_similar_6mers(const size_t, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, float smax0, float smax1, float smax2, float smax3, float smax4, float smax5, int* sort, float* corr);
		void _add_similar_7mers(const size_t, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, int kmer6, float smax0, float smax1, float smax2, float smax3, float smax4, float smax5, float smax6, int* sort, float* corr);

		void _print_kmer(char*);
		void _print_diags(size_t*, const size_t);
		void print_scoring_matrix(size_t imax, size_t jmax);


		//faster methods for generating similar k-mers with presorted list of dimers
		void _add_similar_6mers_fast(const size_t pos, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, float* smax, float* corr);
		void _add_similar_4mers_fast(const size_t pos, int kmer0, int kmer1, int kmer2, int kmer3, float* smax, float* corr);
		
		//based on a 4-mer matrix
		void _add_similar_6mers_from_matrix(const size_t seq_pos, int* kmer, float* corr);
		void _add_similar_4mers_from_matrix(const size_t seq_pos, int* kmer, float* corr);

		inline size_t kmer_idx( const char *kmer ){
			size_t ret   = 0;
			size_t *ptr=aa_powers;
			const char *ch=kmer;
			for(size_t i=0; i<k; ++i){
				ret += (*ptr)*(*ch);
				++ch;
				++ptr;
			}
			return ret;
		}
		
	};
#endif
