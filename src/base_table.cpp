/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "base_table.h"

Base_table::Base_table(	const size_t dbsize, 
			const size_t k, 
			Matrix *m):
			AMINOACID_DIM( m->get_aa_dim() ),
			_int2aa( m->get_int2aa() ),
			_aa2int( m->get_aa2int() ),
			_dbsize(dbsize),
			_k(k){
	_aa_powers = new size_t[k];
	_tablesize = 1;
	for(size_t i=0; i<k; ++i){ 
		_aa_powers[i] = _tablesize;
		_tablesize   *= AMINOACID_DIM;
	}
	const size_t _4mb = 4194304;
	_H                 = new HeapWrapper<_simple_rkmer>( (_4mb/sizeof(_simple_rkmer)) -1 );
	_table             = new _simple_rkmer*[_tablesize];
	_last_in_table     = new _simple_rkmer*[_tablesize];
	_seq_idx2i         = new size_t[_dbsize+1 ];
	_seq_idx2ii        = new size_t[_dbsize+1 ];
	_i2seq_idx         = new size_t[_dbsize+1 ];
	_ii2seq_idx        = new size_t[_dbsize+1 ];
	_max_kmer_score    = new float[_dbsize+1 ];
	_score             = new float[_dbsize+1 ];
	_local_seq_idx2original_seq_idx = new size_t[ dbsize+1 ];
	_lengths           = new float[ dbsize+1 ];
	reset();
}

Base_table::~Base_table(){
	delete [] _aa_powers;
	delete [] _table;
	delete [] _last_in_table;
	delete [] _seq_idx2i;
	delete [] _seq_idx2ii;
	delete [] _i2seq_idx;
	delete [] _ii2seq_idx;
	delete [] _max_kmer_score;
	delete [] _score;
	delete [] _local_seq_idx2original_seq_idx;
	delete [] _lengths;
	delete _H;
}

void Base_table::reset(){
	std::cout << "Base_table::reset\n";
	_set_zero();
	_H->reset();
	_i                     = 1;
	_ii                    = 1;
	_next_free_local_index = 1;
	_saved_rkmers          = 0;
}

size_t Base_table::get_memory_usage(){
	size_t mem=0;
	mem += _tablesize*sizeof(_simple_rkmer*);
	mem += _tablesize*sizeof(_simple_rkmer*);
	mem += 5*(_dbsize+1)*sizeof(size_t);
	mem += 2*(_dbsize+1)*sizeof(float);
	mem += _H->get_instance_memory_usage();
	return mem;
}

size_t Base_table::get_current_memory_usage(){
	size_t mem=0;
	mem += _tablesize*sizeof(_simple_rkmer*);
	mem += _tablesize*sizeof(_simple_rkmer*);
	mem += (5*(_dbsize+1)*sizeof(size_t));
	mem += (2*(_dbsize+1)*sizeof(float));
	if (_H != 0){
		mem += _H->get_current_memory_usage();
	}
	return mem;
}

size_t Base_table::get_number_of_rkmers(){
	return _saved_rkmers;
}

Base_table::MatchIterator Base_table::begin(	const float l, 
						const float min_score,
						const float cov) { 
	return MatchIterator( this, 1, l, min_score, cov );  }
						
Base_table::MatchIterator Base_table::end(){ 
	return MatchIterator( this, _i, 0.0f, 0.0f, 0.0f ); 
}

Base_table::BasicMatchIterator Base_table::begin_iter_all(){ 
	return Base_table::BasicMatchIterator(this, 1);  
}

Base_table::BasicMatchIterator Base_table::end_iter_all() { 
	return Base_table::BasicMatchIterator(this, _i); 
}

void Base_table::add_representative(Sequence *s){
	if(s->length()<_k) return;
	const size_t L             = s->length() - _k;
	const char* const ptr      = s->get_sequence();

	_local_seq_idx2original_seq_idx[_next_free_local_index] =  s->index();
	unsigned int s_index = _next_free_local_index++;
	_lengths[s_index]    = s->length();

	size_t at           = kmer_idx(ptr);
	size_t last         = ptr[0];
	const size_t offset = _k-1;
	const size_t pow    = _aa_powers[offset];
	_simple_rkmer *tmp  = _H->get_next();
	tmp->seq_idx        = s_index;
	tmp->next           = 0;
	++_saved_rkmers;

	if( _table[at]==0 ) _table[at]               = tmp;
	else                _last_in_table[at]->next = tmp;
	_last_in_table[at] = tmp;

	for(size_t j=1; j<=L; ++j){
		at = (at-last)/AMINOACID_DIM;
		at += ptr[j+offset]*pow;	
		last = ptr[j];
		_simple_rkmer *tmp = _H->get_next();
		tmp->seq_idx       = s_index;
		tmp->next          = 0;
		++_saved_rkmers;

		if( _table[at]==0 ) _table[at]               = tmp;
		else                _last_in_table[at]->next = tmp;
		_last_in_table[at] = tmp;
	}
}



void Base_table::add_representative_spaced(Sequence *s){
	// example
	/*                   number of residues in common with first pattern
	XX0XXX0X             -
	 XX0XXX0X            3
	  XX0XXX0X           3
	   XX0XXX0X          3
	    XX0XXX0X         3
	     XX0XXX0X        1
	      XX0XXX0X       1
	       XX0XXX0X      1
	*/

	if(s->length()<8) return;
	const size_t L             = s->length() - 8;
	const char* const ptr      = s->get_sequence();
	_local_seq_idx2original_seq_idx[_next_free_local_index] =  s->index();
	unsigned int s_index = _next_free_local_index++;
	_lengths[s_index]    = s->length();
	char kmer[_k];

	for(size_t j=0; j<=L; ++j){

		kmer[0] = s->get_aa_entry_at(j);
		kmer[1] = s->get_aa_entry_at(j+1);
		kmer[2] = s->get_aa_entry_at(j+3);
		kmer[3] = s->get_aa_entry_at(j+4);
		kmer[4] = s->get_aa_entry_at(j+5);
		kmer[5] = s->get_aa_entry_at(j+7);

		size_t at    = kmer_idx(kmer);
		_simple_rkmer *tmp = _H->get_next();
		tmp->seq_idx       = s_index;
		tmp->next          = 0;
		++_saved_rkmers;
		if( _table[at]==0 ) _table[at]               = tmp;
		else                _last_in_table[at]->next = tmp;
        	_last_in_table[at] = tmp;
	}
}

size_t Base_table::kmer_idx( const char *kmer ){
	size_t ret        = 0;
	size_t *ptr       = _aa_powers;
	const char    *ch = kmer;
	for(size_t i=0; i<_k; ++i){
		ret += (*ptr)*(*ch);
		++ch;
		++ptr;
	}
	return ret;
}

int Base_table::match( Kmer *kmer ){

	for( size_t h=1; h<_i; ++h )_seq_idx2i[_i2seq_idx[h]]=0;
	Kmer_seed_list::Iterator it    = kmer->begin();
	Kmer_seed_list::Iterator end   = kmer->end();
	size_t last_seed_pos = 0;
	_i                   = 1;
	_ii                  = 1;
	number_of_matches    = 0;	
	if( it!=end ) last_seed_pos = (*it)->position; 
	while(it!=end){
		Kmer_seed_list::Kmerseed *seed = (*it);
		if( seed->position!=last_seed_pos){
			for( size_t h=1; h<_ii; ++h ){
				const size_t ii2s_h = _ii2seq_idx[h];
				size_t *s2i_ii2s_h = &_seq_idx2i[ii2s_h];
				if( *s2i_ii2s_h==0 ){
					*s2i_ii2s_h = _i;
					_i2seq_idx[_i]             = ii2s_h;
					_score[_i]                 = _max_kmer_score[h];
					++_i;
				}else{
					_score[ *s2i_ii2s_h ] += _max_kmer_score[h];
				}
				_seq_idx2ii[ii2s_h] = 0;
			}
			last_seed_pos = seed->position;
			_ii = 1;
		}
		_simple_rkmer *rkmer = _table[ seed->index ];
		const float seed_score = seed->score;
		while( rkmer!=0 ){
			++number_of_matches;
			const size_t rkmer_seq_idx = rkmer->seq_idx;
			size_t *rkmer_seq_idx2ii = &_seq_idx2ii[rkmer_seq_idx];
			if(  *rkmer_seq_idx2ii==0 ){
				 *rkmer_seq_idx2ii = _ii;
				_ii2seq_idx[_ii]            = rkmer_seq_idx;
				_max_kmer_score[_ii]        = seed_score;
				++_ii;
			}else{
				if(  seed_score > _max_kmer_score[ *rkmer_seq_idx2ii ] ){
					_max_kmer_score[ *rkmer_seq_idx2ii ] = seed_score;
				}
			}
			rkmer = rkmer->next;
		}
		++it;	
	}
	for( size_t h=1; h<_ii; ++h ){
		if( _seq_idx2i[_ii2seq_idx[h]]==0 ){
			_seq_idx2i[_ii2seq_idx[h]] = _i;
			_i2seq_idx[_i]             = _ii2seq_idx[h];
			_score[_i]                 = _max_kmer_score[h];
			++_i;
		}else{
			_score[_seq_idx2i[_ii2seq_idx[h]]] += _max_kmer_score[h];
		}
		_seq_idx2ii[_ii2seq_idx[h]] = 0;
	}
	return number_of_matches;
}


void Base_table::count_idents( Kmer *kmer ){
	for( size_t h=1; h<_i; ++h ) _seq_idx2i[_i2seq_idx[h]]=0;
	Kmer_seed_list::Iterator it    = kmer->begin();
	Kmer_seed_list::Iterator end   = kmer->end();
	size_t last_seed_pos = 0;
	_i                   = 1;
	if( it!=end ) last_seed_pos = (*it)->position; 
	number_of_matches = 0;
	while(it!=end){
		Kmer_seed_list::Kmerseed *seed = (*it);
		_simple_rkmer *rkmer = _table[ seed->index ];
		while( rkmer!=0 ){
			++number_of_matches;
			const size_t rkmer_seq_idx = rkmer->seq_idx;
			if( _seq_idx2i[rkmer_seq_idx]==0 ){
				_seq_idx2i[rkmer_seq_idx] = _i;
				_i2seq_idx[_i]            = rkmer_seq_idx;
				_score[_i]                = 1;
				++_i;
			}else{
				_score[_seq_idx2i[rkmer_seq_idx]] += 1;
			}
			rkmer = rkmer->next;
		}
		++it;	
	}
}

std::ostream& Base_table::print_kmer(size_t index, std::ostream &out){
	for(size_t i=0; i<_k; ++i){
		size_t aa = index%21;
		out << _int2aa[aa];
		index = index-aa;
		index /= AMINOACID_DIM;
	}
	return out;
}

std::ostream& Base_table::print_table(std::ostream &out){
	for(size_t i=0; i<_tablesize; ++i){
		if( _table[i]!=0 ){
			_simple_rkmer *tmp = _table[i];
			out << std::setw(15) << i << " = ";
			print_kmer(i, out);
			out << " : ";
			while( tmp!=0 ){
				out << "(" << tmp->seq_idx << ") ";
				tmp = tmp->next;
			}
			out << std::endl;
		}
	}
	return out;
}

std::ostream& Base_table::print_matches(std::ostream &out){
	for(size_t i=1; i<_i; ++i){
		out << "       Match with sequence:" << _local_seq_idx2original_seq_idx[ _i2seq_idx[i] ] << " " ;
		out << "Score:" << _score[i] << std::endl;
	}
	return out;
}

void Base_table::write_table_to_file( const char *fn ) throw (std::exception){
	const size_t sof_size_t = sizeof(size_t);
	const size_t sof_uint   = sizeof(unsigned int);
	size_t number_of_used_kmers = 0;
	for( size_t i=0; i<_tablesize; ++i)
		if( _table[i]!=0 ) ++number_of_used_kmers;
	std::ofstream out(fn, std::ios::out | std::ios::binary);
	if(!out) throw MyException("Cannot open '%s' for writing!", fn);
	out.write((char*)&_next_free_local_index, sof_size_t);
	out.write((char*)&number_of_used_kmers, sof_size_t);
	for( size_t i=0; i<_tablesize; ++i){
		if( _table[i]!=0 ){
			size_t count = 0;
			_simple_rkmer *tmp = _table[i];
			while( tmp!=0 ){
				++count;
				tmp = tmp->next;
			}
			out.write((char*)&i, sof_size_t);
			out.write((char*)&count, sof_size_t);
			tmp = _table[i];
			while( tmp!=0 ){
				out.write((char*)&(tmp->seq_idx), sof_uint);
				tmp = tmp->next;
			}
		}
	}
	out.write( (char*) _local_seq_idx2original_seq_idx, sof_size_t*(_next_free_local_index) );
	out.write( (char*) _lengths, sizeof(float)*(_next_free_local_index) );
	out.close();
}


void Base_table::load_from_file( const char *fn ) throw (std::exception){
	reset();
	const size_t sof_size_t = sizeof(size_t);
	const size_t sof_uint   = sizeof(unsigned int);
	std::ifstream in(fn, std::ios::in | std::ios::binary);
	if(!in) throw MyException("Cannot open '%s' for reading!", fn);
	size_t number_of_used_kmers=0;
	in.read((char*)&_next_free_local_index, sof_size_t);
	in.read((char*)&number_of_used_kmers, sof_size_t);
	size_t used_kmers=0;
	while( used_kmers<number_of_used_kmers ){
		size_t i,count;
		in.read((char*)&i, sof_size_t);
		in.read((char*)&count, sof_size_t);
		size_t c=0;
		while( c<count ){
			unsigned int seq_idx;
			in.read((char*)&seq_idx, sof_uint);
			_simple_rkmer *tmp  = _H->get_next();
			tmp->seq_idx = seq_idx;
			tmp->next    = 0;

			if( _table[i]==0 ) _table[i]               = tmp;
			else               _last_in_table[i]->next = tmp;
			_last_in_table[i] = tmp;
			++c;
		}
		++used_kmers;
	}
	in.read( (char*) _local_seq_idx2original_seq_idx, sof_size_t*(_next_free_local_index) );
	in.read( (char*) _lengths, sizeof(float)*(_next_free_local_index) );
	in.close();
}

void Base_table::get_best_match(size_t &index, float &score, float seq_len){
	index = 0;
	score = 0.0f;
	for(size_t i=1; i<_i; ++i){
		if( std::min(seq_len, _lengths[_i2seq_idx[i]]) / std::max(seq_len, _lengths[_i2seq_idx[i]]) < 0.8f ){ 
			continue;
		}
		if( _score[i]/seq_len > score ){
			score = _score[i]/seq_len;
			index = _local_seq_idx2original_seq_idx[ _i2seq_idx[i] ];
		}
	}
}

float Base_table::get_first_match(){
	if( _i==1 ) return 0.0f;
	return _score[1];
}

size_t Base_table::get_kmer_match_count(){ 
	return number_of_matches; 
}

