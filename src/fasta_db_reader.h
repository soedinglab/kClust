#ifndef CM_FASTA_DB_READER_H
#define CM_FASTA_DB_READER_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Reader for protein sequence databases in FastA format

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
//#include <algorithm>
#include <cmath>
#include <climits>
#include <list>

#include "sequence.h"
#include "packed_sequence.h"
#include "my_exception.h"
#include "profile.h"
#include "converter.h"

class Fasta_db_reader{
	public:

		class SeekSortStruct{
			public:
				SeekSortStruct():len(0), pos(0){};
				size_t len;
				size_t index;
				std::streampos pos;
		};

		Fasta_db_reader(std::string fastadb, Matrix *matrix, bool _sc, bool profile_query, std::string inputdir = "", std::string dbsorted_old = "", std::string h_old_file = "");
//		Fasta_db_reader(Params * p, Matrix *matrix);
		~Fasta_db_reader();

		size_t get_index_of_next_sequence();

//		std::ostream& analyze(std::ostream&) throw (std::exception);

		bool move_forward_by(const size_t);
		void move_to(const size_t) throw (std::exception);
		
		Sequence* get_next(bool dummy, bool without_header, bool without_matrix = false);
		Sequence* get( const size_t idx, bool dummy, bool without_header, bool without_matrix = false ) throw (std::exception);
		Profile* get_profile( const size_t pos, bool dummy, bool without_header, bool without_matrix ) throw (std::exception);
		
		Packed_sequence* get_next_packed_without_header();


		Sequence* get_rand();

		bool has_next();
		void reset() throw (std::exception);
		void reset_headers() throw (std::exception);
	
		void seek_sort(const std::string) throw (std::exception);

		std::ostream& write_random_partition(std::ostream&, size_t sequences);

		size_t get_sequence_count(){ return seq_count; }
		
		size_t get_orig_db_size(){ return orig_db_size; }

		std::ostream& analyze_bins(std::ostream &out, int width);

		const size_t get_memory_usage() const;

		const double get_character_count() const{ return chars_d; }
		
		const float get_sum_log_template_len() const{ return sum_log_template_len; }
		
		size_t* get_idx2orig(){return idx2orig;}
		
		size_t* get_orig2idx(){return orig2idx;}
		
		char* get_header_for (int orig_idx);

	private:
	
		Sequence* get_next_sequence (bool dummy, bool without_header);
		Profile* get_next_profile(bool dummy, bool without_header, bool without_matrix);
		
		Sequence* get_sequence( const size_t pos, bool dummy, bool without_header ) throw (std::exception);
		float** get_profile_matrix(int index, int seqlength);
		void _init_profile_settings_();
		void set_profile_members();

		bool _profile_query;
		std::string file;
		std::string alndir;
		std::string clusters_old_file;
		std::string headers_old_file;
		std::ifstream in;
		std::ifstream in_h;
		std::streampos _start_of_data;
		std::streampos _start_of_data_h;
		char *sbuffer; 
		char *hbuffer;
		char *h_old_buffer;
		char *pbuffer;
		bool sc;
		
		const size_t AMINOACID_DIM;
		Matrix * m;
/*		const int*  const aa2int;
		const char* const int2aa;
		const float * const p_back;*/
		//position of the next sequence in the idx2orig array (current_idx = next index for non-profile queries)
		size_t current_idx;
		size_t max_seq_len;
		size_t min_seq_len;
		size_t seq_count;
		size_t orig_db_size;
		// indicates if a sequence with an id i is a member of a profile or not
		size_t *orig2idx;
		// indices of sequences in the order of the database 
		size_t *idx2orig;
		// old clusters represented as a array of lists for efficient access
		// orig_idx -> list of orig_idx
		std::list<size_t> *clusters_old;
		double chars_d;
		float sum_log_template_len;
		std::streampos *_seeks_;
		std::streampos *_hseeks_;	
		SeekSortStruct** _merge_seek_sort(SeekSortStruct*, size_t len);
		void _init_limits_() throw (std::exception);
		// init positions in fasta db file where this sequence begins
		// pos -> seek position
		void _init_seeks_();
		void _init_hseeks_();
};

#endif
