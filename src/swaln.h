#ifndef CM_ALN_H
#define CM_ALN_H

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <vector>

#include "my_exception.h"
#include "sequence.h"
#include "pair_aln.h"

//#include <emmintrin.h>

class SWaln{

	public:

		struct Features{
			float score;
			int len_x;
			int len_y;
			int len_aln;
			int len_equ;
			int idents;
			int pos_matches;
			int neg_zero_matches;
			int gaps;
		};

		struct Alignment{
			char *x;
			char *y;
			Features f;
		};

		enum METHOD{
			OVERLAP_FULL=0, 
			OVERLAP_FLIP=1, 
			OVERLAP_BAND=2, 
			LOCAL_FLIP=3,
			LOCAL_BAND=4
		};	

		//constructor
		SWaln(	METHOD, 
			Sequence*, 
			Sequence*, 
			Matrix *mat,
			const float gapopen=12.0, 
			const float gapextend=2.0, 
			const int diagonal=0, 
			const unsigned int bandwidth=0, 
			const bool v=false);
		//destructor
		~SWaln();	

		//calls an alignment method
		void align();

		//print alignment
		//calls traceback method
		std::ostream& print(std::ostream&, const size_t width=60);

		//compute score, idents,gaps,... only - does not create ASCII strings of the alignment
		void compute_features(Features&);

		//do traceback, save alignment in the reference to the structure
		void traceback(Alignment&);
	
		Pair_aln* traceback_ij();

		void print_traceback(Alignment&);


	private:
		//the alignment method
		const METHOD method;	

		// lengths of sequence x,y 
		const int L_x, L_y;

		// pointer to sequence x,y
		const char* const x;
		const char* const y;

		const int AMINOACID_DIM;
		const char* const int2aa;
		const int* const aa2int;
		const float** const matrix;

		// gapopen and gapextension penalty 
		const float G, E;

		//bandwidth and diagonal
		const int D;
		const int W;
		
		//const float** const matrix;
		//const char* const int2aa;
		
		const bool verbose;

		//score matrix, insertion matrix x, insertion matrix y, traceback matrix
		float **M, **I_x, **I_y;
		char **T; 

		//results of traceback are saved in these vars
		int max_i, max_j, max_j_band;
		float score;		


		//for debugging purposes
		void _init_overlap_matrices_full() throw (std::exception);
		void _align_overlap_full()         throw (std::exception);
		void _delete_overlap_matrices_full();

		void _init_overlap_matrices_flip() throw (std::exception);
		void _align_local_flip()           throw (std::exception);
		void _align_overlap_flip()         throw (std::exception);
		void _delete_overlap_matrices_flip();

		void _init_overlap_matrices_band_flip()             throw (std::exception);
		void _align_local_band_flip(const int, const int)   throw (std::exception);
		void _align_overlap_band_flip(const int, const int) throw (std::exception);
		void _delete_overlap_matrices_band_flip();

		void _traceback_overlap(Alignment&);
		void _traceback_overlap_band(Alignment&);
		void _traceback_local(Alignment&);
		void _traceback_local_band(Alignment&);

		float _compute_overlap_score(char*, char*, size_t);

		void _print_matrix(char**);

		std::ostream& _print_memory_usage(std::ostream&);

		size_t _xstart, _ystart;


};

#endif
