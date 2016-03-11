#ifndef CM_BASE_TABLE_H
#define CM_BASE_TABLE_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//index-table data structure
//implements k-mer scoring match method
//Section 3.1.2.

#include <cstdio>
#include <iostream>
#include <fstream>
#include "matrix.h"
#include "kmer.h"
#include "heap_wrapper.h"

class Base_table{

	public:
		friend class MatchIterator;

		class _simple_rkmer{
			public:
				_simple_rkmer(	const unsigned int seq_idx=0, _simple_rkmer *next=0):
					next(next), seq_idx(seq_idx){}
				_simple_rkmer *next;
				size_t seq_idx;
		};

		class MatchIterator{
			public:
				MatchIterator(	Base_table *bt, 
						const uint h, 
						const float l, 
						const float min_score,
						const float cov):
						_bt(bt),	
						_l(l),
						_min_score(min_score),
						_cov(cov),
						_h(h){
					while( _h<_bt->_i ){
						if( _bt->_score[_h]/_l < _min_score ){
							++_h; 
							continue;
						}
						float tmp = _bt->_lengths[ _bt->_i2seq_idx[_h] ];
						if( tmp>_l ) tmp = _l/tmp;
						else        tmp = tmp/_l;
						if( tmp<_cov ){
							++_h;
							continue;
						}
						break;
					}
				};
				MatchIterator(const MatchIterator &it):
					_bt(it._bt),
					_l(it._l),
					_min_score(it._min_score),
					_cov(it._cov),
					_h(it._h){
				};
				~MatchIterator(){};
				MatchIterator& operator=(const MatchIterator &it){
					this->_h  = it._h;
					this->_bt = it._bt;
					return (*this);
				}
				bool operator==(const MatchIterator &it){
					return this->_h == it._h;
				}
				bool operator!=(const MatchIterator &it){
					return this->_h != it._h;
				}
				MatchIterator& operator++(){
					++_h;
					while( _h<_bt->_i ){
						if( _bt->_score[_h]/_l < _min_score ){
							++_h; 
							continue;
						}
						float tmp = _bt->_lengths[ _bt->_i2seq_idx[_h] ];
						if( tmp>_l ) tmp = _l/tmp;
						else        tmp = tmp/_l;
						if( tmp<_cov ){
							++_h;
							continue;
						}
						break;
					}
					return (*this);
				}
				MatchIterator operator++(int){
					MatchIterator tmp(*this);
					++(*this);
					return tmp;
				}
				size_t get_sequence_index(){
					return _bt->_local_seq_idx2original_seq_idx[_bt->_i2seq_idx[_h]];
				}
				float get_score(){
					return _bt->_score[_h];	
				}
				float get_length(){
					return _bt->_lengths[_h];
				}	
			
			private:
				Base_table *_bt;
				
				const float _l;
				const float _min_score;
				const float _cov;
				uint _h;	
		};

		class BasicMatchIterator{
			public:
				BasicMatchIterator( Base_table *bt, const uint h ):_bt(bt),_h(h){};
				BasicMatchIterator( const BasicMatchIterator &it ):_bt(it._bt),_h(it._h){};
				~BasicMatchIterator(){};
				BasicMatchIterator& operator=(const BasicMatchIterator &it){
					this->_h  = it._h;
					this->_bt = it._bt;
					return (*this);
				}
				bool operator==(const BasicMatchIterator &it){ return this->_h == it._h; }
				bool operator!=(const BasicMatchIterator &it){ return this->_h != it._h; }
				BasicMatchIterator& operator++(){
					++_h;
					return (*this);
				}
				BasicMatchIterator operator++(int){
					BasicMatchIterator tmp(*this);
					++(*this);
					return tmp;
				}
				size_t get_sequence_index(){ 
					return _bt->_local_seq_idx2original_seq_idx[_bt->_i2seq_idx[_h]]; 
				}
				float get_score(){ return _bt->_score[_h]; }
				float get_length(){ return _bt->_lengths[_h]; }	
			private:
				Base_table *_bt;
				uint _h;
					
		};

		MatchIterator begin(	const float l, 
					const float min_score,
					const float cov);
		MatchIterator end();

		BasicMatchIterator begin_iter_all();
		BasicMatchIterator end_iter_all();

		Base_table(const size_t dbsize, const size_t k, Matrix *m);
		~Base_table();

		void add_representative( Sequence* );
		void add_representative_spaced( Sequence* );

		int match( Kmer* );
		size_t get_kmer_match_count();

		void count_idents( Kmer* );

		void get_best_match(size_t &index, float &score, float seq_len);	

		float get_first_match();
		
		// set all indexes to zero according to k-mer similarity scoring algorithm (algorithm 3 in section 3.1.2.)
		void reset();
		size_t kmer_idx( const char *kmer );
		
		void write_table_to_file(const char *fn ) throw (std::exception);
		void load_from_file(const char *fn ) throw (std::exception);

		size_t get_memory_usage();
		size_t get_number_of_rkmers();
		size_t get_current_memory_usage();

		std::ostream& print_kmer(size_t, std::ostream&);
		std::ostream& print_table(std::ostream&);
		std::ostream& print_matches(std::ostream&);
	private:
		const size_t      AMINOACID_DIM;
		const char* const _int2aa;
		const int*  const _aa2int;
		const size_t      _dbsize;
		const size_t      _k;
		
		size_t           _tablesize;
		size_t          *_aa_powers;
		
		size_t           _saved_rkmers;
		//array of pointers to simple-lists - initially zero
		_simple_rkmer **_table;
		_simple_rkmer **_last_in_table;

		size_t *_seq_idx2i;
		size_t *_i2seq_idx;
		size_t _i;
	
		size_t *_seq_idx2ii;
		size_t *_ii2seq_idx;
		size_t _ii;

		float *_max_kmer_score;
		float *_score;

		size_t number_of_matches;

		//each representative is mapped to local index
		//the range is 1..(number of sequences in the table) 
		size_t *_local_seq_idx2original_seq_idx;

		//counter holds the number of the next free local-index for a representative
		//counter-1 is the current number of representatives in the table
		size_t _next_free_local_index;

		float *_lengths;
		HeapWrapper<_simple_rkmer> *_H;

		inline void _set_zero(){
			_simple_rkmer **ptr = _table;
			_simple_rkmer **end = _table + _tablesize;
			while( ptr!=end){ *ptr++=0; }
			for(size_t i=0; i<=_dbsize; ++i){
				_seq_idx2i[i]  = 0;
				_seq_idx2ii[i] = 0;
				_local_seq_idx2original_seq_idx[i] = 0;
				_lengths[i] = 0.0f;
			}
		}
};

#endif
