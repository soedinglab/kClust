#ifndef CM_PACKED_SEQUENCE_H 
#define CM_PACKED_SEQUENCE_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Implements a sequence container that stores amino acid sequences in a compressed form (6aa in 4bytes)
//used for merging large sequence databases and filtering of identities using a hash 

#include <inttypes.h>
#include <iostream>
#include "sequence.h"
#include "my_exception.h"

class Packed_sequence{
	public:
		Packed_sequence();
		Packed_sequence( 	const char* const header_src, 
								const char* const sequence_src, 
								const size_t index_src=0 );
		Packed_sequence( const Sequence * const) throw (std::exception);
		Packed_sequence( const Packed_sequence &other );
		~Packed_sequence();
		Packed_sequence& operator=(const Packed_sequence &other);
		bool operator==( const Packed_sequence &other );
		std::ostream& print( std::ostream &out);
		inline const unsigned char * const get_seq_ptr() const { return (unsigned char*)sequence; }
		inline const size_t length() const { return l;    }
		inline const size_t byte_length() const { return l*sizeof(uint32_t);    }
		inline const size_t aa_length() const { return al; }
		inline const uint32_t get_checksum(){ return checksum; }
		inline const char* const get_header() const{ return header; }
		inline const size_t index() const{ return idx; }
		const size_t get_memory_usage() const;
	private:
		void init( 	const char* const header_src, 
						const char* const sequence_src, 
						const size_t index );
		const uint32_t compute_checksum() const;
		uint32_t *sequence;
		char *header;
		size_t l;
		size_t al;
		size_t idx;
		uint32_t checksum;
};

#endif
