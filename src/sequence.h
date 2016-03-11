#ifndef CM_SEQUENCE_H
#define CM_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//A container class for protein sequences

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>

#include "my_exception.h"
#include "matrix.h"

//4MB
#define READ_HBUFFER_SIZE 4194304
//1MB
#define READ_SBUFFER_SIZE 1048576
//20MB
#define READ_PBUFFER_SIZE 20971520

class Matrix;
class Sequence{
	
	public:
		enum format{
			fasta     = 0,
			clustal   = 1, //not implemented
			stockholm = 2 //not implemented
		};

	public:

		Sequence(	const char* const header, 
				const char* const sequence, 
				bool _sc,
				size_t index=0,
				Matrix *m=NULL,
				bool dummy=true
			) throw (std::exception);

/*		Sequence(	const char *filename, 
				const format&, 
				Matrix *m, 
				size_t index=0
			) throw (std::exception);*/

		//create random seq of length 'length' with respect to background amino acid frequencies
		Sequence(	const size_t length,
				Matrix *m, 
				size_t index=0
			) throw (std::exception);


		virtual ~Sequence();

		//inline functions HAVE TO reside in header files to be accessible from other cpp files
		inline const char * const get_sequence() const { return _seq; }
		inline const char * const get_header()   const { return _header; }
		inline char *sequence()                        { return _seq; }
		inline const size_t length() const             { return len; }
		// sequence identifier
		inline const size_t index()  const             { return idx; }
		inline const bool is_dummy() const             { return primitive; }
		void set_index(const size_t);

		std::ostream& print_debug(std::ostream&);
		std::ostream& write(std::ostream&, const size_t width=60);
		std::ostream& write_with_pseudo_header(std::ostream&, const size_t width=60);
		std::ostream& write_sequence(std::ostream&, const size_t width=60);
		// clears profile data in Profile objects
		virtual void clear(){};
		virtual size_t get_memory_usage();

		static size_t get_length( const char* const ptr, const size_t l);
		void replace_J();
		void rescale_matrix() throw (std::exception);
		Matrix * const get_matrix() { return m;}
		
		virtual const float** const get_scoring_matrix();
		virtual const float * const get_max_scores();
		virtual inline void get_kmer_at (int * kmer, int k, int pos) { for (int i = 0; i < k; i++) *(kmer++) = _seq[pos++]; }
		virtual inline int get_kmer_entry_at (int pos) {  return _seq[pos];}
		virtual inline int get_aa_entry_at (int pos) { return _seq[pos]; }
		virtual float get_score_with (int seq_pos, char aa);
		virtual float* get_score_correction();
		
	protected:
		virtual void calc_score_correction (int pseudocounts = 20);
/*		const size_t AMINOACID_DIM;
		const int*  const aa2int;
		const char* const int2aa;*/
		Matrix * m;
		// primitive: aa characters in _seq
		// !primitive: aa numbers in _seq
		bool primitive;
		float* correction;
		float* zero_correction;

		char *_header;
		// (char) aa2int of the original amino acids in the sequence
		char *_seq;
		size_t len;
		size_t _hlen;
		size_t idx;
		bool sc;
		
	private:
		void readFasta(const char*) throw (std::exception);
		void _init(	const char* const header, 
						const size_t hlen, 
						const char* const sequence, 
						const size_t slen) throw (std::exception);
		
};


#endif

