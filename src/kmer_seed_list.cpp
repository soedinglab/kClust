/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "kmer_seed_list.h"

Kmer_seed_list::Kmer_seed_list(){
	H = new Kmer_seed_list::KmerseedHeap( __MEMB_SIZE );
	reset();
}

Kmer_seed_list::~Kmer_seed_list(){
	delete H;
}

Kmer_seed_list::Iterator::Iterator(	Kmerseed *seed=0, 
												size_t load=0, 
												size_t current=0):
												seed(seed), 
												load(load), 
												current(current){}

Kmer_seed_list::Iterator::~Iterator(){}

size_t Kmer_seed_list::get_instance_memory_usage(){
	//std::cerr << "SEED BLOCKS:" << H->get_block_count() << std::endl;
	return H->get_instance_memory_usage() + 20;
}

void Kmer_seed_list::append(const uint index, const float score, const ushort position){
	Kmerseed *tmp = H->get_next();
	tmp->index    = index;
	tmp->score    = score;
	tmp->position = position;
	++_load;
}

void Kmer_seed_list::reset(){
	_load = 0;
	H->reset();
}

Kmer_seed_list::Iterator Kmer_seed_list::begin(){ 
	return Kmer_seed_list::Iterator(H->get_first(), _load, 0);
}

Kmer_seed_list::Iterator Kmer_seed_list::end(){ 
	return Kmer_seed_list::Iterator(0, _load, _load);
}

size_t Kmer_seed_list::load(){return _load;}

Kmer_seed_list::Iterator& Kmer_seed_list::Iterator::operator=(const Iterator& other){
	seed    = other.seed;
	load    = other.load;
	current = other.current;
	return (*this);
}

bool Kmer_seed_list::Iterator::operator==(const Kmer_seed_list::Iterator& other){ 
	return (current == other.current); 
}

bool Kmer_seed_list::Iterator::operator!=(const Kmer_seed_list::Iterator& other){ 
	return (current != other.current); 
}

Kmer_seed_list::Iterator& Kmer_seed_list::Iterator::operator++(){
	seed = seed->next;
	++current;	
	return(*this);
}

Kmer_seed_list::Iterator Kmer_seed_list::Iterator::operator++(int){
	Iterator tmp(*this);
	++(*this);
	return(tmp) ;
}

Kmer_seed_list::Kmerseed* Kmer_seed_list::Iterator::operator*(){ 
	return seed; 
}

