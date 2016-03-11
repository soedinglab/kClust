/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "recycle_table.h"

Recycle_table::Recycle_table(const uint dbsize, const uint k, Matrix *m):
	AMINOACID_DIM( m->get_aa_dim() ),
	_int2aa( m->get_int2aa() ),
	_aa2int( m->get_aa2int() ),
	_k(k), 
	_dbsize(dbsize){

	_aa_powers = new uint[k];

	_tablesize = 1;
	for(size_t i=0; i<k; ++i){ 
		_aa_powers[i] = _tablesize;
		_tablesize   *= AMINOACID_DIM;
	}

	_table           = new _rkmer*[_tablesize];
	_last_in_table   = new _rkmer*[_tablesize];

	_used_table_entries = new size_t[_tablesize];
	_table_occ = 0;

	const size_t _4mb = 4194304;
	Rkmer_heap       = new HeapWrapper<_rkmer>( (_4mb/sizeof(_rkmer)) -1 );
	Match_heap       = new HeapWrapper<_match>( (_4mb/sizeof(_match)) -1);

	_seq_idx2i       = new uint[dbsize+1];
	_seq_idx2ii      = new uint[dbsize+1];
	_i2seq_idx       = new uint[dbsize+1];
	_ii2seq_idx      = new uint[dbsize+1];

	_max_kmer_score  = new float[dbsize+1];
	_max_kmer_pos    = new int[dbsize+1];
	_score           = new float[dbsize+1];

	_kmer_match      = new _match*[dbsize+1];
	_last_kmer_match = new _match*[dbsize+1];

	_local_seq_idx2original_seq_idx = new uint[dbsize+1];

	_set_zero();
	_i                     = 1;
	_ii                    = 1;
	_next_free_local_index = 1;
}	

Recycle_table::~Recycle_table(){
	delete Rkmer_heap;
	delete Match_heap;
	delete [] _table;
	delete [] _last_in_table;
	delete [] _used_table_entries;
	delete [] _aa_powers;
	delete [] _seq_idx2i;
	delete [] _seq_idx2ii;
	delete [] _i2seq_idx;
	delete [] _ii2seq_idx;
	delete [] _max_kmer_score;
	delete [] _max_kmer_pos;
	delete [] _score;
	delete [] _kmer_match;
	delete [] _last_kmer_match;
	delete [] _local_seq_idx2original_seq_idx;
}

size_t Recycle_table::get_memory_usage(){
	size_t mem=0;
	mem += _tablesize*sizeof(_rkmer*);
	mem += _tablesize*sizeof(_rkmer*);
	mem += _tablesize*sizeof(size_t);
	mem += 6*(_dbsize+1)*sizeof(size_t);
	mem += 2*(_dbsize+1)*sizeof(float);
	mem += 2*(_dbsize+1)*sizeof(_match*);
	mem += Rkmer_heap->get_instance_memory_usage();
	mem += Match_heap->get_instance_memory_usage();
	return mem;
}

uint Recycle_table::get_number_of_matching_seqs(){ return _i-1; }

void Recycle_table::reset(){
	for( size_t h=0; h<_table_occ; ++h) _table[_used_table_entries[h]] = 0; 
	_table_occ = 0;
	for( size_t h=1; h<_i; ++h ) _seq_idx2i[_i2seq_idx[h]]=0;
	Rkmer_heap->reset();
	Match_heap->reset();
	_i                     = 1;
	_ii                    = 1;
	_next_free_local_index = 1;
}

void Recycle_table::add_representative( Sequence *s ){
	if(s->length()<_k) return;
	const uint L       = s->length()-_k;
	const char *ptr    = s->get_sequence();
	_local_seq_idx2original_seq_idx[_next_free_local_index] =  s->index();
	const uint seq_idx = _next_free_local_index++;
	const size_t pow   = _aa_powers[_k-1];
	uint offset        = _k-1;
	uint last          = ptr[0];
	uint at            = kmer_idx(ptr);
	_rkmer *tmp        = Rkmer_heap->get_next();
	tmp->seq_idx       = seq_idx;
	tmp->pos           = 0;
	tmp->next          = 0;

	if( _table[at]==0 ){ 
		_table[at]               = tmp;
		_used_table_entries[_table_occ++] = at;
	}else                _last_in_table[at]->next = tmp;
        _last_in_table[at] = tmp;

	for(size_t j=1; j<=L; ++j){
		at   = (at-last)/AMINOACID_DIM;
		at  += ptr[j+offset]*pow;	
		last = ptr[j];

		tmp = Rkmer_heap->get_next();
		tmp->seq_idx = seq_idx;
		tmp->pos     = j;
		tmp->next    = 0;

		if( _table[at]==0 ){
			_table[at]               = tmp;
			_used_table_entries[_table_occ++] = at;
		}else                _last_in_table[at]->next = tmp;
        	_last_in_table[at] = tmp;
	}
}

void Recycle_table::add_representative_spaced( Sequence *s ){
	if(s->length()<_k) return;
	const uint L       = s->length()-_k-1;
	const char *ptr    = s->get_sequence();
	_local_seq_idx2original_seq_idx[_next_free_local_index] =  s->index();
	const uint seq_idx = _next_free_local_index++;
	char kmer[_k];
	for(size_t j=0; j<=L; ++j){
//		memcpy(kmer, &ptr[j], 2*sizeof(char));
//		memcpy(kmer+2, &ptr[j]+3, (_k-2)*sizeof(char));

		kmer[0] = s->get_aa_entry_at(j);
		kmer[1] = s->get_aa_entry_at(j+1);
		kmer[2] = s->get_aa_entry_at(j+3);
		kmer[3] = s->get_aa_entry_at(j+4);

		size_t at    = kmer_idx(kmer);
		_rkmer *tmp  = Rkmer_heap->get_next();
		tmp->seq_idx = seq_idx;
		tmp->pos     = j;
		tmp->next    = 0;

		if( _table[at]==0 ){ 
			_table[at]               = tmp;
			_used_table_entries[_table_occ++] = at;
		}else   _last_in_table[at]->next = tmp;
		
		_last_in_table[at] = tmp;
	}
}

std::ostream& Recycle_table::print_table(std::ostream &out){	
	for(uint i=0; i<_tablesize; ++i){
		if( _table[i]!=0 ){
			_rkmer *tmp = _table[i];
			out << i << " = ";
			print_kmer(i, out);
			out << " : ";
			while( tmp!=0 ){
				out << "(" << tmp->seq_idx << "," << tmp->pos << ") ";
				tmp = tmp->next;
			}
			out << std::endl;
		}
	}
	return out;
}

std::ostream& Recycle_table::print_kmer(size_t index, std::ostream &out){
	for(size_t i=0; i<_k; ++i){
		size_t aa = index%21;
		out << _int2aa[aa];
		index = index-aa;
		index /= AMINOACID_DIM;
	}
	return out;
}

std::ostream& Recycle_table::print_matches(std::ostream &out){
	for(size_t i=1; i<_i; ++i){
		out << "Match with sequence:" << _i2seq_idx[i] << " ";
		out << "Score:" << _score[i] << std::endl;
	}
	return out;
}

int Recycle_table::match(Kmer *kmer) throw (std::exception) {
	for( size_t h=1; h<_i; ++h ) _seq_idx2i[_i2seq_idx[h]]=0;
	// iterators over similar k-mers for the given k-mer
	Kmer_seed_list::Iterator it  = kmer->begin();
	Kmer_seed_list::Iterator end = kmer->end();
	size_t last_seed_pos         = 0;
	_i                           = 1;
	_ii                          = 1;
	if( it!=end ) last_seed_pos = (*it)->position; 
	// seed: sim. k-mer index, score with the query k-mer, position in the query sequence
	Kmer_seed_list::Kmerseed *seed;
	_rkmer *rkmer;
	_match *match;
	int db_match_num = 0;
	while(it!=end){
		seed = (*it);
		if( seed->position!=last_seed_pos){
			for( size_t h=1; h<_ii; ++h ){
				if( _seq_idx2i[_ii2seq_idx[h]]==0 ){
					_seq_idx2i[_ii2seq_idx[h]]  = _i;
					_i2seq_idx[_i]              = _ii2seq_idx[h];
					_score[_i]                  = _max_kmer_score[h];
					match                       = Match_heap->get_next();
					match->i                    = last_seed_pos;
					match->j                    = _max_kmer_pos[h];
					match->score                = _max_kmer_score[h];
					match->next                 = 0;
					_kmer_match[_i]             = match;
					_last_kmer_match[_i]        = match;
					++_i;
				}else{
					const size_t iii            = _seq_idx2i[_ii2seq_idx[h]];
					_score[iii]                += _max_kmer_score[h];
					match                       = Match_heap->get_next();
					match->i                    = last_seed_pos;
					match->j                    = _max_kmer_pos[h];
					match->score                = _max_kmer_score[h];
					match->next                 = 0;
					_last_kmer_match[iii]->next = match;
					_last_kmer_match[iii]       = match;
				}
				_seq_idx2ii[_ii2seq_idx[h]] = 0;
			}
			last_seed_pos = seed->position;
			_ii = 1;
		}
		
		rkmer = _table[seed->index];
		while( rkmer!=0 ){
			// first match for the sequence seq_idx
			if( _seq_idx2ii[rkmer->seq_idx]==0 ){
				_seq_idx2ii[rkmer->seq_idx] = _ii;
				_ii2seq_idx[_ii]            = rkmer->seq_idx;
				_max_kmer_score[_ii]        = seed->score;
				_max_kmer_pos[_ii]          = rkmer->pos;
				++_ii;
			}else{
				if( seed->score > _max_kmer_score[_seq_idx2ii[rkmer->seq_idx]] ){
					_max_kmer_score[_seq_idx2ii[rkmer->seq_idx]] = seed->score;
					_max_kmer_pos[_seq_idx2ii[rkmer->seq_idx]]   = rkmer->pos;
				}
			}
			rkmer = rkmer->next;
			db_match_num++;
		}
		++it;	
	}

	for( size_t h=1; h<_ii; ++h ){
		if( _seq_idx2i[_ii2seq_idx[h]]==0 ){
			_seq_idx2i[_ii2seq_idx[h]]  = _i;
			_i2seq_idx[_i]              = _ii2seq_idx[h];
			_score[_i]                  = _max_kmer_score[h];
			match                       = Match_heap->get_next();
			match->i                    = last_seed_pos;
			match->j                    = _max_kmer_pos[h];
			match->score                = _max_kmer_score[h];
			match->next                 = 0;
			_kmer_match[_i]             = match;
			_last_kmer_match[_i]        = match;

			++_i;
			
		}else{
			const size_t iii            = _seq_idx2i[_ii2seq_idx[h]];
			_score[iii]                += _max_kmer_score[h];
			match                       = Match_heap->get_next();
			match->i                    = last_seed_pos;
			match->j                    = _max_kmer_pos[h];
			match->score                = _max_kmer_score[h];
			match->next                 = 0;
			_last_kmer_match[iii]->next = match;
			_last_kmer_match[iii]       = match;
		}
		_seq_idx2ii[_ii2seq_idx[h]] = 0;
	}
	return db_match_num;
}

int Recycle_table::match_full(Kmer *kmer) throw (std::exception){

	for( size_t h=1; h<_i; ++h )
		_seq_idx2i[_i2seq_idx[h]]=0;

	Kmer_seed_list::Iterator it    = kmer->begin();
	Kmer_seed_list::Iterator end   = kmer->end();
	_i                             = 1;
	Match_heap->reset();

	Kmer_seed_list::Kmerseed *seed;
	_rkmer *rkmer;
	_match *match;
	int db_match_num = 0;
	while(it!=end){
		seed = (*it);
		rkmer = _table[seed->index];
		while( rkmer!=0 ){
			if( _seq_idx2i[rkmer->seq_idx]==0 ){
				_seq_idx2i[rkmer->seq_idx]  = _i;
				_i2seq_idx[_i]              = rkmer->seq_idx;
				match                       = Match_heap->get_next();
				match->i                    = seed->position;
				match->j                    = rkmer->pos;
				match->score                = seed->score;
				match->next                 = 0;
				_kmer_match[_i]             = match;
				_last_kmer_match[_i]        = match;
				++_i;
			}else{
				const size_t iii            = _seq_idx2i[rkmer->seq_idx];
				match                       = Match_heap->get_next();
				match->i                    = seed->position;
				match->j                    = rkmer->pos;
				match->score                = seed->score;
				match->next                 = 0;
				_last_kmer_match[iii]->next = match;
				_last_kmer_match[iii]       = match;
			}
			rkmer = rkmer->next;
			db_match_num++;
		}
		++it;	
	}
	return db_match_num;
}

Recycle_table::MatchIterator Recycle_table::begin() { return MatchIterator(this, 1);  }
Recycle_table::MatchIterator Recycle_table::end()   { return MatchIterator(this, _i); }

Recycle_table::_match *Recycle_table::get_best_match() throw (std::exception) {
	float maxs  = 0.0f;
	size_t maxp = 0;
	for(size_t h=1; h<_i; ++h){
		if( _score[h]>maxs ){
			maxs = _score[h];
			maxp = h;
		}
	}
	if( _i==1 ){
		throw MyException("No kmer matches -> lower threshold for similar k-mers!");
	}
	return _kmer_match[maxp];
}
