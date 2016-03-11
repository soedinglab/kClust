#ifndef CM_ALIGNMENT_H
#define CM_ALIGNMENT_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION
//Container class for multiple sequence alignments
//Implements sequence weights and consensus sequence computation

#include <string>
#include <vector>
#include <utility>
#include "matrix.h"
#include "simple_hash.h"

#define CM_ALIGNMENT_BUFFER_SIZE 4194304

class Alignment{

	public:
		enum format{
			fasta     = 0,
			clustal   = 1,
			a3m       = 2
		};
		Alignment( const char *filename, const format fileformat, Matrix *m ) throw (std::exception);
		~Alignment();
		void init_max_entropy_weights(); //Henikoff&Henikoff '94
		std::ostream& write(	std::ostream &out, 
									const format fileformat, 
									const size_t seq_width=80, 
									const size_t header_width=30);
		std::string get_consensus_sequence( int line_width=-1, int method=1 );	

	private:
		Matrix *matrix;
		const int AMINOACID_DIM;
		const int*  const aa2int;
		const char* const int2aa;
		const std::string filename;
		size_t columns;
		size_t number_of_sequences;
		std::string *headers;
		char **sequences;
		float *weights;
		int *residue_counts;
		void _read_fasta_alignment( const char *filename ) throw (std::exception);
		void _process_data( 	const char hbuffer[], 
					const char sbuffer[], 
					std::vector<std::string> &headers_vector, 
					std::vector<char*> &sequences_vector,
					std::vector< std::pair<int,int> > &sequences_length_vector ) throw (std::exception);
		std::ostream& _write_clustal(std::ostream &out, const size_t seq_width); 
		std::ostream& _write_fasta(std::ostream &out, const size_t seq_width);
		std::ostream& _write_a3m(std::ostream &out, const size_t seq_width);
		
		inline float fast_log2(float x);

};

#endif
