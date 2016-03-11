/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "packed_sequence.h"

Packed_sequence::Packed_sequence():sequence(0), header(0), l(0), al(0), checksum(0){
}

Packed_sequence::Packed_sequence( const Sequence * const s ) throw (std::exception){
	if( !s->is_dummy() ) 
		throw MyException("Invalid sequence object provided for packed_sequence constructor!");
	init( s->get_header(), s->get_sequence(), s->index() );
}

Packed_sequence::Packed_sequence( 	const char* const header_src, 
												const char* const sequence_src, 
												const size_t index_src ){
	init( header_src, sequence_src, index_src );
}

void Packed_sequence::init( 	const char* const header_src, 
										const char* const sequence_src,  
										const size_t index_src ){
	idx = index_src;
	const size_t hlen = strlen(header_src)+1;
	header = new char[hlen];
	memcpy( header, header_src, hlen*sizeof(char) ); 
	const size_t full_len = strlen(sequence_src);
	const size_t aalen = Sequence::get_length( sequence_src, full_len );
	sequence = new uint32_t[ (aalen/6)+1 ];
	char to_int[256];
	for( int i='A'; i<'A'+26; ++i ) to_int[i] = i-'A';
	uint32_t powers[6];
	powers[0] = 1;
	for( int i=1; i<6; ++i ) powers[i] = powers[i-1]*26;
	uint32_t _6mer_index = 0;
	l                    = 0;
	int _6mer_pos        = 0;
	for( size_t i=0; i<full_len; ++i ){
		if( !isalpha( sequence_src[i] ) ) continue;
		//std::cerr << powers[_6mer_pos]*(sequence_src[i]-'A') << std::endl;
		_6mer_index += powers[_6mer_pos]*(sequence_src[i]-'A');
		++_6mer_pos;
		if( _6mer_pos==6 ){
			sequence[l++] = _6mer_index;
			_6mer_index   = 0;
			_6mer_pos     = 0;
		}
	}
	al = 6*l + _6mer_pos;	
	if( _6mer_pos!=0 ) sequence[l++] = _6mer_index;
	checksum = 0;
}

Packed_sequence::Packed_sequence( const Packed_sequence &other ){
	l = other.l;
	al = other.al;
	sequence = new uint32_t[l];
	memcpy( sequence, other.sequence, sizeof( uint32_t )*l );
	const size_t hlen = strlen(other.get_header())+1;
	header = new char[hlen];
	memcpy( header, other.get_header(), hlen*sizeof(char) ); 
}	

Packed_sequence::~Packed_sequence(){
	delete [] header;
	delete [] sequence;
}

Packed_sequence& Packed_sequence::operator=(const Packed_sequence &other){
	delete [] header;
	delete [] sequence;
	this->l        = other.l;
	this->al       = other.al;
	this->checksum = other.checksum;
	const size_t hlen = strlen(other.get_header())+1;
	header = new char[hlen];
	memcpy( header, other.get_header(), hlen*sizeof(char) ); 

	sequence = new uint32_t[l];
	memcpy( sequence, other.sequence, sizeof( uint32_t )*l );

	return (*this);
}

bool Packed_sequence::operator==( const Packed_sequence &other ){
	if( other.al!=al ) return false;
	for( size_t i=0; i<l; ++i ) 
		if( other.sequence[i]!=sequence[i] ) 
			return false;
	return true;	
}

std::ostream& Packed_sequence::print( std::ostream &out ){
	out << header << std::endl;
	uint32_t powers[6];
	powers[0] = 1;
	for( int i=1; i<6; ++i ) powers[i] = powers[i-1]*26;
	size_t chars=0;
	for( size_t i=0; i<l; ++i ){
			uint32_t _6mer_index = sequence[i];
			for( int j=0; j<6; ++j ){
				char c = _6mer_index%26;
				out << char('A'+c);
				_6mer_index -= c;
				_6mer_index /= 26;
				++chars;
				if( chars==al) break;
			}
	}
	out << std::endl;
	return out;
}

const uint32_t Packed_sequence::compute_checksum() const {
	uint32_t ret=0;
	for( size_t i=0; i<l; ++i ) ret ^= sequence[i];
	return ret;
}

const size_t Packed_sequence::get_memory_usage() const{
	size_t ret=0;
	ret += (sizeof(size_t) + ((strlen(header)+1)*sizeof(char))); 
	ret += (sizeof(size_t) + l*sizeof(uint32_t)); 
	ret += sizeof(Packed_sequence);
	return ret;
}
