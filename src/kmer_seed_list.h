#ifndef CM_KMER_SEED_LIST_H
#define CM_KMER_SEED_LIST_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

#include <iostream>
#include "heap_wrapper.h"
#define __MEMB_SIZE 10000

class Kmer_seed_list{

	public:
		typedef unsigned int uint;
		typedef unsigned short ushort;

		class Kmerseed{
			public:
				Kmerseed *next;
				uint index;
				float score;
				ushort position;

				Kmerseed(){
					index    = 0;
					score    = 0.0f;
					position = 0;
					next     = 0;
				}
		
				Kmerseed(	const uint index,
								const float score,
								const ushort position,
								Kmerseed *next=0	):
								next(next), 
								index(index), 
								score(score), 
								position(position){}

				~Kmerseed(){}
		};

		class KmerseedHeap : public HeapWrapper<Kmerseed>{
			public:
				KmerseedHeap(const size_t blocksize):HeapWrapper<Kmerseed>(blocksize){
					_link_block( _first->data );	
				}
				virtual ~KmerseedHeap(){
					//DESTRUCTOR OF PARENT CLASS IS AUTOMATICALLY CALLED
				}								
				inline Kmerseed* get_next(){
					if( _i==_BLOCK_SIZE ){
						_block *tmp = _current->next;
						if( tmp==0){
							tmp            = new _block(_BLOCK_SIZE, _bc++);
							tmp->next      = 0;
							_current->next = tmp;
							_link_block( tmp->data );
							_current->data[_BLOCK_SIZE-1].next = tmp->data;
						}
						_current = tmp;
						_i       = 0;
					}
					return &_current->data[_i++];
				}
				Kmerseed* get_first(){ return _first->data; }
				void reset(){ _i=0; _current=_first; }
			private:
				inline void _link_block(Kmerseed *data){
					Kmerseed *last      = data;
					Kmerseed *it        = &data[1];
					const Kmerseed *end = data + _BLOCK_SIZE;
					while( it!=end){
						last->next = it;
						last       = it++;
					}
					last->next=0;
				}
		};

		class Iterator{
			public:
				friend class Kmerseed;
				Iterator(Kmerseed*, size_t, size_t);
				~Iterator();
				Iterator& operator=(const Iterator&);
				bool operator==(const Iterator&);
				bool operator!=(const Iterator&);
				Iterator& operator++();
				Iterator operator++(int);
				Kmerseed* operator*();
			private:
				Kmerseed *seed;
				size_t load, current;
		};

		Kmer_seed_list();
		~Kmer_seed_list();
		void append(const uint index, const float score, const ushort position);
		void reset();
		Iterator begin();
		Iterator end();
		size_t load();
		size_t get_instance_memory_usage();
		
	private:
		//number of used elements
		unsigned int _load;
		unsigned int _index_of_next_free;
		KmerseedHeap *H;
};

#endif
