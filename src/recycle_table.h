#ifndef CM_RECYCLE_TABLE_H 
#define CM_RECYCLE_TABLE_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//index-table data structure
//implements k-mer matching method for computation of match lists 
//stores position i in seq. x and j in  seq. y for each match

#include "kmer.h"
#include "heap_wrapper.h"

class Recycle_table{

	public:
		typedef unsigned int uint;	

		friend class MatchIterator;

		// k-mers object, stores seq-id and the position
		class _rkmer{
			public:
				_rkmer(){}
				_rkmer(	const uint seq_idx, 
					const uint position,
					_rkmer *next=0,
					_rkmer *global_next=0):
					seq_idx(seq_idx), 
					pos(position),
					next(next){}
				uint seq_idx;
				uint pos;
				_rkmer *next;
		};
		
		//match-list element
		class _match{
			public:
				_match(){};
				unsigned short i;
				unsigned short j;
				float score;
				_match *next;	
		};

		//Iterator for match-lists
		class MatchIterator{
			public:
				MatchIterator(Recycle_table *rt=0, uint h=0){
					this->h  = h;
					this->rt = rt;
				};
				MatchIterator(const MatchIterator &it){
					this->h  = it.h;
					this->rt = it.rt;
				};
				~MatchIterator(){};
				MatchIterator& operator=(const MatchIterator &it){
					this->h  = it.h;
					this->rt = it.rt;
					return (*this);
				}
				bool operator==(const MatchIterator &it){
					return this->h == it.h;
				}
				bool operator!=(const MatchIterator &it){
					return this->h != it.h;
				}
				MatchIterator& operator++(){
					//if( h<=rt->_i ) 
					++h;
					return (*this);
				}
				MatchIterator operator++(int){
					MatchIterator tmp(*this);
					++(*this);
					return tmp;
				}
				uint get_sequence_index(){
					return rt->_local_seq_idx2original_seq_idx[rt->_i2seq_idx[h]];
				}
				float get_score(){
					return rt->_score[h];	
				}				
				_match *get_match_list(){
					return rt->_kmer_match[h];
				}			
			private:
				uint h;
				Recycle_table *rt;
		};

		Recycle_table(uint dbsize, uint k, Matrix*);
		~Recycle_table();

		//reset the table - allocated memory is reused
		void reset();

		//add the k-mers of a representative
		void add_representative(Sequence*);

		//add spaced k-mers of a representative
		//PATTERN: XX0X...
		void add_representative_spaced(Sequence*);

		inline uint kmer_idx( const char *kmer ){
			size_t ret        = 0;
			unsigned int *ptr = _aa_powers;
			const char *ch    = kmer;
			for(size_t i=0; i<_k; ++i){
				ret += (*ptr)*(*ch);
				++ch;
				++ptr;
			}
			return ret;
		}	

		std::ostream& print_kmer(size_t, std::ostream&);
		std::ostream& print_table(std::ostream&);
		std::ostream& print_matches(std::ostream &out);

		//create match-lists of maximum scoring k-mers 
		int match(Kmer*) throw (std::exception);

		//create match-lists of ALL k-mers 
		int match_full(Kmer*) throw (std::exception);

		//this method works only with a preceeding call to match()
		//returns a pointer to the list of the maximum scoring representative
		_match* get_best_match() throw (std::exception);

		uint get_number_of_matching_seqs();

		MatchIterator begin();
		MatchIterator end();

		size_t get_memory_usage();

	private:
		const size_t AMINOACID_DIM;
		const char* const _int2aa;
		const int*  const _aa2int;

		const uint _k, _dbsize;
		
		_rkmer **_table, **_last_in_table;

		size_t *_used_table_entries;
		size_t _table_occ;

		uint _tablesize;
		uint *_aa_powers;
		
		uint *_seq_idx2i;
		uint *_i2seq_idx;
		uint _i;
	
		uint *_seq_idx2ii;
		uint *_ii2seq_idx;
		uint _ii;

		//each representative gets a local index
		uint *_local_seq_idx2original_seq_idx;

		uint _next_free_local_index;

		float *_max_kmer_score;
		int   *_max_kmer_pos;
		float *_score;
		
		//arrays with pointers to lists in the index-table
		//_kmer_match points to the first element of each MATCH-LIST of the query and a representative
		_match **_kmer_match;
		//_last_kmer_match points to the last element of each MATCH-LIST of the query and a representative
		_match **_last_kmer_match;

		//allocators for list objects
		HeapWrapper<_rkmer> *Rkmer_heap;
		HeapWrapper<_match> *Match_heap;

		inline void _set_zero(){
			_rkmer **ptr = _table;
			_rkmer **end = _table + _tablesize;
			while( ptr!=end){ *ptr++=0; }

			for(uint i=0; i<=_dbsize; ++i){
				_seq_idx2i[i]  = 0;
				_seq_idx2ii[i] = 0;
			}
		}
};




#endif
