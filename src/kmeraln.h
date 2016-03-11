#ifndef CM_KMERALN_H
#define CM_KMERALN_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION
//Implements k-mer dynamic programming and the k-mer graph data structure

#include <vector>
#include <cmath>
#include "kmer.h"
#include "recycle_table.h"
#include "matrix.h"
#include "pair_aln.h"
#include "heap_wrapper.h"

//#define __CM_DEBUG_KMERALN 1
//#define __CM_DEBUG_KMERALN_MATCH_LIST 1
//#define __CM_DEBUG_KMERALN_DYNAMIC 1
//#define __CM_DEBUG_KMERALN_TRACEBACK 1

class Kmeraln{
	typedef unsigned int uint;

	public:
		enum Method{
			FULL_BAND=0,
			FAST_ADDR=1
		};

		class _extended_match{
			public:
			int i;
			int j;
			
			float kmerscore;
			bool dummy;
			float dynscore;
			_extended_match *inext;
			//_extended_match *gnext;
			_extended_match *intra_d_next;
			_extended_match *intra_d_prev;
			_extended_match *inter_d_next;
			_extended_match *inter_d_prev;
			_extended_match *traceback;
			_extended_match(int i=0, 
					int j=0, 
					float score=0.0f,
					bool dummy=false):
					i(i), 
					j(j), 
					kmerscore(score),
					dummy(dummy){}
			~_extended_match(){}
		};

		Kmeraln(	const int k, 
				Matrix *mat,
				const float p_m,
				const int delta=50,
				const float G=12.0f, 
				const float E=2.0f, 
				const float F=0.1f,
				const float Lamda=0.3f, 
				const float K=0.09f, 
				const float H=0.34f,
				const double db_len=1,
				const float sum_log_db_len=0, 
				const float db_seqs=1,
				const int   KMER_OFFSET=-1,
				const float COV=-INFINITY ) throw (std::exception);

		~Kmeraln();

		size_t get_memory_usage();

		// m     Specifies the algorithm for dynamic programming 
 		// mlist List of kmer matches, i: position in sequence x, j: position in sequence y
		// x     Sequence for which the kmer-seeds were generated.
		// y     Sequence for which 'identity' kmers  are saved in the matching-table
		void align(	Method, 
				Recycle_table::_match *mlist, 
				Sequence *x, 
				Sequence *y) throw (std::exception);
		
		std::ostream& print(
					std::ostream&,
					Sequence *x, 
					Sequence *y, 
					const float** const matrix, 
					const size_t width=60 ) throw (std::exception);

		std::ostream& print_real_trace( std::ostream&) throw (std::exception);
		
		std::ostream& print_alignment(std::ostream& out) throw (std::exception);

		float get_score();
		float get_rel_score();
		float get_evalue();
		float get_aln_cov_short();
		float get_aln_cov_long();

		size_t get_number_of_matches();

		std::ostream& print_trace(std::ostream&);
		std::ostream& print_vertex_list(std::ostream&);
		
		void set_coverage_criterion(float coverage);

		Pair_aln *get_aln();
		
	private:
		
		const int AMINOACID_DIM;
		const char* const int2aa;

		const float p_m;
		const int delta;

		//flag for artificial matches
		bool _add_grid_kmers_flag;	

		// gap penalties	
		const float G, E, F;
		
		// statistical parameters
		const float Lamda, K, H;
	
		// Database statistic
		const double db_len;
		const float sum_log_db_len, db_seqs;
		
		// drop-off score
		const float _dropoff;

		//grid kmer match score
		const float _grid_score;

		//offset per kmer in traceback
		//set to 1 for 4mers and 3mers
		//set to 2 for spaced 4mers
		int _kmer_offset;

		//kmer length	
		const int _k;

		//coverage	
		float COV;

		//sequence lengths
		int Lx, Ly;

		//number of fast addresses - sizeof(diag_address[])
		int C;
		int W;

		//pointers to sequences
		char *X;
		char *Y;
		Sequence *X_seq;
		Sequence *Y_seq;

		//score of kmeraln and coverage
		float _score;
		float aln_cov_short;
		float aln_cov_long;

		//vertex with maximum k-mer aln score
		Kmeraln::_extended_match *_trace;

		//arrays for optimization of gaps -> traceback
		float *forward, *backward;

		//fast addressing array
		_extended_match **diag_ptr;

		//pointer to the first vertex in the graph
		_extended_match *_first;

		//pointer to alignment
		Pair_aln *aln;

		//memory allocator for vertices
		HeapWrapper<_extended_match> *ExtMatchHeap;

		void _init_fast_addr() throw (std::exception);
		void _delete_fast_addr();

		//copy list of kmer-matches into the vertices of the graph
		void _init_list(Recycle_table::_match*) throw (std::exception);

		//k-mer dynamic programming with fast addressing
		void _fill_fast();
		
		//k-mer dynamic programming without fast addressing
		//starts either at the bottom or top of the Delta band to find possible predecessors
		void _fill();
		
		//check if 'candidate' is in the Delta window 
		inline bool _in_window(_extended_match *current, _extended_match *candidate);

		//optimize k-mer alignment
		Pair_aln* traceback();

		//optimize gap
		inline float _compute_gapx(	const int, 
						const int, 
						const int, 
						const int,
						Pair_aln *alnp, 
						int &alnpos);

		inline float _compute_gapy(	const int, 
						const int, 
						const int, 
						const int,
						Pair_aln *alnp, 
						int &alnpos);

		//computes the score of an alignment
//		float _compute_overlap_score(char*, char*, size_t, const float**);

};

#endif
