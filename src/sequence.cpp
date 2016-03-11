/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "sequence.h"


Sequence::Sequence(	const char* const header, 
			const char* const sequence, 
			bool _sc,
			size_t index,
			Matrix *matrix,
			bool dummy) throw (std::exception){
	correction = NULL;
	zero_correction = NULL;
	sc = _sc;
	if (!dummy){
		m = matrix;
		primitive = false;
		len = 0;
		idx = index;
		const size_t hlen = strlen(header);
		const size_t slen = strlen(sequence);
		_init(header, hlen, sequence, slen);
	}
	else{ //dummy
		m = matrix;
		primitive = true;
		len =0;
		idx = index;
		const size_t hlen = strlen(header);
		const size_t slen = strlen(sequence);

		//allocate memory for header- and sequence-data 
		_header = new char[hlen+1];
		_seq    = new char[slen+1];

		// copy data + appending '\0' into buffers
		memcpy(_header, header, hlen+1);
		memcpy(_seq, sequence, slen+1);

		//determine number of amino acids in the sequence
		len = Sequence::get_length( _seq, slen );
	}
}

/*Sequence::Sequence(	const char* fn, 
			const Sequence::format &f, 
			Matrix *matrix,
			size_t index) throw (std::exception):
			m(matrix),
			primitive(false),
			len(0),
			idx(index){
	switch(f){
		case Sequence::fasta : readFasta(fn); break;
		default              : throw MyException("Unknown sequence format!\n");
	}
}*/

Sequence::Sequence( 	const size_t length,
			Matrix *matrix,  
			const size_t index) throw (std::exception):
/*			AMINOACID_DIM(m->get_aa_dim()),
			aa2int(m->get_aa2int()), 
			int2aa(m->get_int2aa()),*/
			m(matrix),
			primitive(false),
			len(length),
			idx(index){
	_header = new char[60];
	_seq    = new char[len];
	sprintf(_header, ">random, length: %zu", len);

	float aa_freq[20];
	float sum=0.0;
	for( size_t i=0; i<20; ++i ){
		sum += m->get_p_back()[i];
		aa_freq[i] = sum;
	}

	for( size_t i=0; i<len; ++i ){
		int aa;
		float r = rand()/(float)RAND_MAX;
		for( size_t j=0; j<20; ++j ){
			aa=j;
			if( r<aa_freq[j] ) break;	
		}
		_seq[i] = aa;
	}
}

void Sequence::set_index(const size_t index){ idx = index; }

size_t Sequence::get_memory_usage(){
	
	size_t ret = 7*sizeof(size_t);
	//header
	ret += sizeof(char)*_hlen;
	//sequence
	ret += sizeof(char)*len;
	return ret;
}

void Sequence::readFasta(const char *fn) throw (std::exception) {
	char *hbuffer = new char[READ_HBUFFER_SIZE];
	char *sbuffer = new char[READ_SBUFFER_SIZE];
	std::ifstream in(fn);
	if( in.fail() ) throw MyException("Cannot read '%s'", fn);
	// move to first header line 
	while( !in.eof() ){
		char c = in.get();
		if( c=='>' ){ in.unget(); break; }	
		if( isspace(c) ) continue;
		throw MyException("File '%s' is not in correct FastA format!", fn);
	}
	//read header and sequence
	if(!in.eof()){	
		//extracts with appended '\0'
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		size_t hb = in.gcount();
		if( hb==(READ_HBUFFER_SIZE-1) ) 
			throw MyException("Error reading '%s'\nLong header line found, increase READ_HBUFFER_SIZE!", fn);	
		if( hbuffer[0] == '>' ){
			in.get(sbuffer, READ_SBUFFER_SIZE, '>');
			size_t sb = in.gcount();
			if( sb==(READ_SBUFFER_SIZE-1) ) 
				throw MyException("Error reading '%s'\nLong sequence found, increase READ_SBUFFER_SIZE!", fn);
			if( sb==0 ) 
				throw MyException("File '%s' is not in correct FastA format!", fn);
			_init(hbuffer, hb, sbuffer, sb);
		}else throw MyException("File '%s' is not in correct FastA format!", fn);
	}else         throw MyException("File '%s' is not in correct FastA format!", fn);

	delete [] hbuffer;
	delete [] sbuffer;
}


Sequence::~Sequence(){
	delete [] _header;
	delete [] _seq;
	if (correction != NULL) delete[] correction;
	if (zero_correction != NULL) delete[] zero_correction;
}

void Sequence::_init(const char* const header, const size_t hlen, const char* const sequence, const size_t slen) throw (std::exception){
	
	const int*  const aa2int = m->get_aa2int();
	// copy header to new memory location so that it takes exactly as much space as it needs
	_header = new char[hlen+1];
	memcpy(_header, header, hlen+1);
	_hlen = hlen;

	// len is the real length of sequence
	len = Sequence::get_length(sequence, slen);
	if(len==0)throw MyException("No sequence data found for '%s'!\n", _header);
	_seq = new char[len];
	
	// convert all characters where required to other characters (for invalid AAs like B) and then to corresponding ints, store them in _seq
	char *ptr = _seq;
	for(size_t i=0; i<slen; ++i){
		if( isalpha(sequence[i]) ){ 
			const char cur = toupper(sequence[i]);
			switch(cur){
				case 'A': *ptr = aa2int['A']; break;
				case 'B': *ptr = aa2int['N']; break;
				case 'C': *ptr = aa2int['C']; break;
				case 'D': *ptr = aa2int['D']; break;
				case 'E': *ptr = aa2int['E']; break;
				case 'F': *ptr = aa2int['F']; break;
				case 'G': *ptr = aa2int['G']; break;
				case 'H': *ptr = aa2int['H']; break;
				case 'I': *ptr = aa2int['I']; break;
				case 'J': *ptr = aa2int['L']; break;
				case 'K': *ptr = aa2int['K']; break;
				case 'L': *ptr = aa2int['L']; break;
				case 'M': *ptr = aa2int['M']; break;
				case 'N': *ptr = aa2int['N']; break;
				case 'O': *ptr = aa2int['X']; break;
				case 'P': *ptr = aa2int['P']; break;
				case 'Q': *ptr = aa2int['Q']; break;
				case 'R': *ptr = aa2int['R']; break;
				case 'S': *ptr = aa2int['S']; break;
				case 'T': *ptr = aa2int['T']; break;
				case 'U': *ptr = aa2int['C']; break;
				case 'V': *ptr = aa2int['V']; break;
				case 'W': *ptr = aa2int['W']; break;
				case 'X': *ptr = aa2int['X']; break;
				case 'Y': *ptr = aa2int['Y']; break;
				case 'Z': *ptr = aa2int['E']; break;
				
				default: 
					std::cerr << "Warning: Ignoring invalid character '";
					std::cerr << sequence[i] << "' at position " << i;
					std::cerr << " in sequence " << header << std::endl;
					continue;
			}
			++ptr;
		}else if( !isspace(sequence[i]) ){
			std::cerr << "Warning: Ignoring invalid character '";
			std::cerr << sequence[i] << "' at position " << i;
			std::cerr << " in sequence " << header << std::endl;
		}
	}
	calc_score_correction();
}

size_t Sequence::get_length( const char* const ptr, const size_t l){
	size_t ret = 0;
	for( size_t i=0; i<l; ++i )if( isalpha(ptr[i]) ) ++ret;
	return ret;
}

std::ostream& Sequence::print_debug(std::ostream &out){
	out << "INDEX     : " << idx       << std::endl;
	out << "LENGTH    : " << len       << std::endl;
	out << "DUMMY     : " << primitive << std::endl;
	out << "HEADER    : " << _header   << "$" << std::endl;
	if(!primitive){
		out << "NUMERIC   : ";
		for(size_t i=0; i<len; ++i) out << (int)_seq[i] << " ";
		out << "$" << std::endl;
	}
	if( !primitive ){
		const char* const int2aa = m->get_int2aa();
		out << "SEQUENCE  : ";
		for(size_t i=0; i<len; ++i) out << int2aa[(int)_seq[i]] << " ";
		out << "$" << std::endl;
	}else out << _seq << "$" << std::endl;
	return out;
}

std::ostream& Sequence::write(std::ostream &out, const size_t width){
	out << _header << std::endl;
	if( !primitive ){
		const char* const int2aa = m->get_int2aa();
		for(size_t i=0; i<len; i+=width){
			for(size_t j=i; j<std::min(i+width, len); ++j)out << int2aa[(int)_seq[j]];	
			out << std::endl;
		}	
	}
	else out << _seq;
	return out;
}

std::ostream& Sequence::write_with_pseudo_header(std::ostream &out, const size_t width){
	out << ">"<< idx << std::endl;
	if( !primitive ){
		const char* const int2aa = m->get_int2aa();
		for(size_t i=0; i<len; i+=width){
			for(size_t j=i; j<std::min(i+width, len); ++j) out << int2aa[(int)_seq[j]];	
			out << std::endl;
		}	
	}
	else out << _seq;
	return out;
}

std::ostream& Sequence::write_sequence(std::ostream &out, const size_t width){
	if( !primitive ){
		const char* const int2aa = m->get_int2aa();
		for(size_t i=0; i<len; i+=width){
			for(size_t j=i; j<std::min(i+width, len); ++j) out << int2aa[(int)_seq[j]];	
			out << std::endl;
		}	
	}
	else out << _seq;
	return out;
}

void Sequence::replace_J(){
	if(primitive){
		const size_t slen = strlen(_seq);
		for( size_t i=0; i<slen; ++i){ 
			if( _seq[i]=='J' )      _seq[i] = 'L';
			else if( _seq[i]=='j' ) _seq[i] = 'l';
			else if( _seq[i]=='O' ) _seq[i] = 'X';
			else if( _seq[i]=='o' ) _seq[i] = 'x';
		}
	}
}

void Sequence::rescale_matrix() throw (std::exception){
	m->wn_rescale(this);
}


const float ** const Sequence::get_scoring_matrix() { 
	return m->get_matrix();
}

const float * const Sequence::get_max_scores() {
	return m->get_diagonal_scores();
}

float Sequence::get_score_with (int seq_pos, char aa) { 
	float* corr;
	if (sc) corr = correction;
	else corr = zero_correction;
	if (!primitive) return (m->get_matrix()[_seq[seq_pos]][aa] + corr[seq_pos]);
	else {
		const int*  const aa2int = m->get_aa2int();
		return (m->get_matrix()[aa2int[_seq[seq_pos]]][aa] + corr[seq_pos]);
	}
}

void Sequence::calc_score_correction(int pseudocounts){
	if (primitive) return;
	
	const int aa_dim = m->get_aa_dim();
	const float* f = m->get_p_back();
	float bg_freq [aa_dim];
	int counts[aa_dim];
	for (int i = 0; i < aa_dim; i++) counts[i] = 0;
	for (int i = 0; i < len; i++) counts[_seq[i]]++;
	for (int i = 0; i < aa_dim; i++){
		bg_freq[i] = (f[i]*((float) pseudocounts) + counts[i]) / (pseudocounts + len);
	}
	
	correction = new float [len];
	const float* m_bg = m->get_p_back();
	for (int i = 0; i < len; i++){
		if (bg_freq[_seq[i]] != 0.0){
			correction[i] = 0.5 * m->_log2(m_bg[_seq[i]] / bg_freq[_seq[i]]);
		}
		else{
			correction[i] = 0.0;
		}
	}
	
	zero_correction = new float [len];
	memset(zero_correction, 0.0f, len*sizeof(float));
}

float* Sequence::get_score_correction(){
	if (sc) return correction;
	return zero_correction;
}
