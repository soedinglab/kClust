/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "block_manager.h"

BlockManager::BlockManager( Fasta_db_reader *db_reader, Clusters *clusters ):
	dbsize(db_reader->get_sequence_count()), 
	db_reader(db_reader), 
	clusters(clusters){
	seqs = new Sequence*[dbsize+1];
	Sequence **ptr = seqs;
	const Sequence **end = (const Sequence**) seqs + dbsize + 1;
	while( ptr!=end ) *ptr++=0;
	current_table_block.set_index( 1 );
	_current_sequence_memory_usage = 0;
}

BlockManager::~BlockManager(){
	free_sequences();
	delete [] seqs;
}

void BlockManager::free_sequences(){
	for( size_t i=0; i<dbsize+1; ++i ){
		if( seqs[i]!=0 ){ 
			delete seqs[i];
			seqs[i] = 0;
		}
	}
	_current_sequence_memory_usage = 0;
}

BlockManager::Iterator BlockManager::begin(){
	return block_list.begin();
}

BlockManager::Iterator BlockManager::end(){
	return block_list.end();
}

BlockManager::Block BlockManager::get_current_table_block(){
	return current_table_block;
}

size_t BlockManager::get_current_table_index(){
	return current_table_block.get_index();
}

void BlockManager::add_block( const Block &b ){
	block_list.append( b );
	current_table_block = b;
}

void BlockManager::save_representative( Sequence *r ){
    r->clear();
    _current_sequence_memory_usage += r->get_memory_usage();
    seqs[r->index()] = r;
}

size_t BlockManager::get_memory_usage(){
	return dbsize*sizeof(Sequence*)+_current_sequence_memory_usage;
}

Sequence* BlockManager::get_representative( size_t idx ){
	return seqs[idx];
}

void BlockManager::load( const Block &b ){
    std::cout << "Loading block, begin: " << b.begin << ", end: " << b.end << "\n";
    free_sequences();
	db_reader->reset();
	db_reader->move_to( b.begin );
	const size_t last = b.end-1;
	while( db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, true, true);
		if( clusters->is_representative_in( s->index(), b.get_index() ) ){
			seqs[ s->index() ] = s;
			if( s->index()==last) break;
		}else{
			if( s->index()==last){
				delete s;
				break;
			}
			delete s;
		}
	}
	if (b.end > current_table_block.end)
		current_table_block = b;
//	std::cout << "Overall memory usage of sequences in block manager after load: " << _current_sequence_memory_usage << "\n\n";
}

size_t BlockManager::get_number_of_blocks(){
	return block_list.length();	
}


