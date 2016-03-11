/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "my_exception.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "alignment.h"

Alignment::Alignment( 	const char *filename, 
			const Alignment::format fileformat, 
			Matrix *m ) throw (std::exception):
			matrix(m),
			AMINOACID_DIM(m->get_aa_dim()),
			aa2int(m->get_aa2int()), 
			int2aa(m->get_int2aa()),
			filename( filename ){
	weights             = 0;
	residue_counts      = 0;
	sequences           = 0;
	headers             = 0;
	number_of_sequences = 0;
	columns             = 0;
	switch(fileformat){
		case(Alignment::fasta): _read_fasta_alignment( filename );
					break;
		default: throw MyException("Unknown file format: '%i'", fileformat);
	}
}

Alignment::~Alignment(){
	for( size_t i=0; i<number_of_sequences; ++i)
		delete [] sequences[i];
	delete [] sequences;
	delete [] headers;
	delete [] weights;
	delete [] residue_counts;
}

void Alignment::_read_fasta_alignment( const char *filename ) throw (std::exception){
	std::ifstream in(filename);
	if( in.fail() ) throw MyException("Cannot read '%s'", filename);
	
	char *hbuffer = new char[CM_ALIGNMENT_BUFFER_SIZE];
	char *sbuffer = new char[CM_ALIGNMENT_BUFFER_SIZE];
	
	// move to first header
	while( !in.eof() ){
		char c = in.get();
		if( c=='>' ){ in.unget(); break; }
		if (c == '#') { in.getline(hbuffer, CM_ALIGNMENT_BUFFER_SIZE, '\n'); continue; }	
		if( isspace(c) ) continue;
		throw MyException("File '%s' is not in correct FastA format!", filename);
	}

	std::vector<std::string> headers_vec;
	std::vector< char* > sequences_vec;
	std::vector< std::pair<int,int> >sequences_length_vec;

	//read header and sequence
	if(!in.eof()){	
		while(!in.eof()){	
			//extracts with appended '\0'
			in.getline(hbuffer, CM_ALIGNMENT_BUFFER_SIZE, '\n');	
			size_t hb = in.gcount();
			if( hb==(CM_ALIGNMENT_BUFFER_SIZE-1) ) 
				throw MyException("Error reading '%s'\nLong header line found!", filename);
			if( hbuffer[0] == '>' ){
				in.get(sbuffer, CM_ALIGNMENT_BUFFER_SIZE, '>');
				size_t sb = in.gcount();
				if( sb==(CM_ALIGNMENT_BUFFER_SIZE-1) ) 
					throw MyException("Error reading '%s'\nLong sequence found!", filename);
				if( sb==0 ) throw MyException("File '%s' is not in correct FastA format!", filename);
				_process_data( hbuffer, sbuffer, headers_vec, sequences_vec, sequences_length_vec );
			}else throw MyException("File '%s' is not in correct FastA format!", filename);
		}
	}else throw MyException("File '%s' is empty!", filename);
	
	if( headers_vec.size()!=sequences_vec.size() )  
		throw MyException("File '%s' is not in correct FastA format!", filename);
	number_of_sequences = sequences_vec.size();

	//check length
	int slen = 0;
	for( size_t i=0; i<number_of_sequences; ++i){
		if( slen<1 ) slen = sequences_length_vec[i].first;
		if( slen==0 ) 
			throw MyException("No sequence available for '%s' (sequence %i)", headers_vec[i].c_str(), i+1);
		if( slen!=sequences_length_vec[i].first ) 
			throw MyException("Sequence %i has different length '%s'", i+1, headers_vec[i].c_str());
		if( sequences_length_vec[i].second==0 ) 
			throw MyException("Sequence %i consists of gaps only! '%s'", i+1, headers_vec[i].c_str());
	}

	columns        = slen;
	headers        = new std::string[number_of_sequences];
	sequences      = new char*[number_of_sequences];
	residue_counts = new int[number_of_sequences];

	for( size_t i=0; i<number_of_sequences; ++i){
		headers[i]        = headers_vec[i];
		sequences[i]      = sequences_vec[i];
		residue_counts[i] = sequences_length_vec[i].second;
	}
	delete [] sbuffer;
	delete [] hbuffer;
}

void Alignment::_process_data( 	const char hbuffer[], 
				const char sbuffer[], 
				std::vector<std::string> &headers_vector, 
				std::vector<char*> &sequences_vector,
				std::vector< std::pair<int,int> > &sequences_length_vector) throw (std::exception){
	headers_vector.push_back( std::string( hbuffer ) );
	const size_t L  = strlen( sbuffer );
	char *tmp       = new char[L];
	size_t pos      = 0;
	size_t residues = 0;
	for( size_t i=0; i<L; ++i ){
		if( isalpha(sbuffer[i]) ){ 
			const char cur = toupper(sbuffer[i]);
			switch(cur){
				case 'A': tmp[pos] = aa2int['A']; break;
				case 'B': tmp[pos] = aa2int['N']; break;
				case 'C': tmp[pos] = aa2int['C']; break;
				case 'D': tmp[pos] = aa2int['D']; break;
				case 'E': tmp[pos] = aa2int['E']; break;
				case 'F': tmp[pos] = aa2int['F']; break;
				case 'G': tmp[pos] = aa2int['G']; break;
				case 'H': tmp[pos] = aa2int['H']; break;
				case 'I': tmp[pos] = aa2int['I']; break;
				case 'J': tmp[pos] = aa2int['L']; break;
				case 'K': tmp[pos] = aa2int['K']; break;
				case 'L': tmp[pos] = aa2int['L']; break;
				case 'M': tmp[pos] = aa2int['M']; break;
				case 'N': tmp[pos] = aa2int['N']; break;
				case 'O': tmp[pos] = aa2int['X']; break;
				case 'P': tmp[pos] = aa2int['P']; break;
				case 'Q': tmp[pos] = aa2int['Q']; break;
				case 'R': tmp[pos] = aa2int['R']; break;
				case 'S': tmp[pos] = aa2int['S']; break;
				case 'T': tmp[pos] = aa2int['T']; break;
				case 'U': tmp[pos] = aa2int['C']; break;
				case 'V': tmp[pos] = aa2int['V']; break;
				case 'W': tmp[pos] = aa2int['W']; break;
				case 'X': tmp[pos] = aa2int['X']; break;
				case 'Y': tmp[pos] = aa2int['Y']; break;
				case 'Z': tmp[pos] = aa2int['E']; break;
				default: 
					throw MyException("Invalid character found: '%c' at position %i in sequence %i\n %s\n", sbuffer[i], i, headers_vector.size(), hbuffer);
					break;
			}
			++pos;
			++residues;
		}else if( sbuffer[i]=='-' ){
			tmp[pos] = AMINOACID_DIM;
			++pos;
		}else if( !isspace(sbuffer[i]) ){
			throw MyException("Invalid character found: '%c' at position %i in sequence %i\n %s\n", sbuffer[i], i, headers_vector.size(), hbuffer);
		}
	}
	char *raw = new char[pos];
	memcpy(raw, tmp, pos);
	delete [] tmp;
	sequences_vector.push_back( raw );
	sequences_length_vector.push_back( std::pair<int,int>(pos, residues) );
}

std::ostream& Alignment::write( std::ostream &out, const format fileformat, const size_t seq_width, const size_t header_width){
	switch(fileformat){
		case( Alignment::a3m )     : return _write_a3m( out, seq_width );
		case( Alignment::fasta   ) : return _write_fasta( out, seq_width );
		case( Alignment::clustal ) : return _write_clustal( out, seq_width);
		default: throw MyException("Unknown file format: '%i'", fileformat);
	}
	return out;
}

std::ostream& Alignment::_write_fasta( std::ostream &out, const size_t seq_width ){
	for( size_t i=0; i<number_of_sequences; ++i){
		out << headers[i] << std::endl;
		for( size_t j=0; j<columns; ++j ){
			if( sequences[i][j]==AMINOACID_DIM ) out << '-';
			else                                 out << int2aa[sequences[i][j]];
			if( (j+1)%seq_width==0 ) out << std::endl;
		}
		if( columns%seq_width!=0 ) out << std::endl;
	}
	return out;
}

std::ostream& Alignment::_write_a3m( std::ostream &out, const size_t seq_width ){
	if( weights==0 ) init_max_entropy_weights();

	float *column_residues_weight = new float[columns];
	for( size_t i=0; i<columns; ++i ) column_residues_weight[i] = 0.0f;

	for( size_t i=0; i<number_of_sequences; ++i)
		for( size_t j=0; j<columns; ++j )
			if( sequences[i][j]<(AMINOACID_DIM-1) ) column_residues_weight[j] += weights[i];

	char buffer[seq_width];
	for( size_t i=0; i<number_of_sequences; ++i){
		out << headers[i] << std::endl;
		size_t printed=0;
		for( size_t j=0; j<columns; ++j ){
			if( sequences[i][j]==AMINOACID_DIM ){ 
				if( column_residues_weight[j]>=0.5f ) buffer[printed++] = '-';
			}else{ 
				if( column_residues_weight[j]>=0.5f ) buffer[printed++] = int2aa[sequences[i][j]];
				else                                  buffer[printed++] = tolower( int2aa[sequences[i][j]] );
			}
			if( printed==seq_width ){ 
				out.write(buffer, seq_width); 
				out << std::endl;
				printed=0;
			}
		}
		if( printed!=seq_width ){
			out.write(buffer, printed);
			out << std::endl;
		}
	}
	delete [] column_residues_weight;
	return out;
}

std::ostream& Alignment::_write_clustal( std::ostream &out, const size_t seq_width ){
	Simple_hash<int> check_headers(number_of_sequences);
	bool shrink_headers=false;
	for( size_t i=0; i<number_of_sequences; ++i ){
		std::string n = headers[i].substr(1, 30);
		if( !check_headers.exists( n.c_str(), 30 ) ){
			check_headers.put( n.c_str(), 30, 1);
		}else{
			shrink_headers = true;
			break;
		}	
	}
	int nr_len = (int)(log10(number_of_sequences))+2; // length of number + '_'
	char hbuffer[32];
	std::stringstream format;
	format << "%-s_%-" << nr_len << "i ";

	std::vector<std::string> clustal_headers;
	for( size_t i=0; i<number_of_sequences; ++i ){
		if( shrink_headers )  sprintf( hbuffer, format.str().c_str(), headers[i].substr(1, 30-nr_len).c_str(), i );
		else	              sprintf( hbuffer, "%-30s ", headers[i].substr(1, 30).c_str() );
		bool trailings=true;
		for( int h=29; h>=0; --h ){
			if( isspace(hbuffer[h]) && !trailings ) hbuffer[h]='_';
			if( !isspace(hbuffer[h]) && trailings ) trailings=false;
		}
		clustal_headers.push_back( std::string(hbuffer) );
	}
	out << "CLUSTAL X (0.00) multiple sequence alignment" << std::endl;
	size_t j=0;
	char buffer[seq_width];
	while( j<columns ){	
		for( size_t i=0; i<number_of_sequences; ++i){
			out << clustal_headers[i];	
			size_t b=0;
			for(size_t jj=j; jj<columns; ++jj ){
				if( sequences[i][jj]<AMINOACID_DIM )
					buffer[b++] = int2aa[sequences[i][jj]];
				else
					buffer[b++] = '-';
				if( b==seq_width ) break;
			}
			out.write(buffer, b);
			out << std::endl;
		}
		j+=seq_width;
		out << std::endl;
	}
	return out;
}

void Alignment::init_max_entropy_weights(){
	weights = new float[number_of_sequences];
	//alloc mem
	const size_t aadim_gap = AMINOACID_DIM+1;
	int **k_x = new int*[columns];
	int *m    = new int[columns];
	for( size_t i=0; i<columns; ++i){
		k_x[i] = new int[aadim_gap];
		m[i]   = 0;
		for( size_t j=0; j<aadim_gap; ++j ) k_x[i][j]=0;
	}
	//count aa frequency per column
	for( size_t i=0; i<number_of_sequences; ++i )
		for( size_t j=0; j<columns; ++j )
			++k_x[ j ][ sequences[i][j] ];
	//init mi[] - # different aa in column i
	const int real_aadim = AMINOACID_DIM-1;
	for( size_t i=0; i<columns; ++i ){
		for( int j=0; j<real_aadim; ++j )
			if( k_x[i][j] ) ++m[i];
		if( m[i]==0 ) m[i]=1;
	}
	//compute weights
	for( size_t k=0; k<number_of_sequences; ++k ){
		weights[k] = 0.0f;
		for( size_t i=0; i<columns; ++i){
			//do not count for X'x
			if( sequences[k][i]<(AMINOACID_DIM-1) )
				weights[k] += 1.0f/( m[i] * k_x[i][ sequences[k][i] ] * residue_counts[k]);
		}
	}
	//normalize
	float sum=0.0f;
	for( size_t k=0; k<number_of_sequences; ++k ) sum += weights[k];
	//std::cerr << "Sum of weights:" << sum << std::endl;
	for( size_t k=0; k<number_of_sequences; ++k ) weights[k]/=sum;

	//free mem
	for( size_t i=0; i<columns; ++i){
		delete [] k_x[i];
	}
	delete [] k_x;
	delete [] m;
}

std::string Alignment::get_consensus_sequence( int line_width, int method ){
	if( weights==0 ) init_max_entropy_weights();
	std::stringstream str;
	//alloc mem
	const int aadim_gap = AMINOACID_DIM+1;
	float **col_freqs = new float*[columns];
	for( size_t i=0; i<columns; ++i){
		col_freqs[i] = new float[aadim_gap];
		for( int j=0; j<aadim_gap; ++j ) col_freqs[i][j] = 0.0f;
	}
	for( size_t k=0; k<number_of_sequences; ++k )
		for( size_t i=0; i<columns; ++i)
		  if( sequences[k][i]!=AMINOACID_DIM ) // throw out by incrising dim of col_freqs[]
				col_freqs[i][ sequences[k][i] ] +=  weights[k];
	const int real_aadim = AMINOACID_DIM-1;
	switch (method) {
		
		case 2:		// max. profile score with consensus
		{
			const float * const freq = matrix->get_p_back();
			for( size_t i=0; i<columns; ++i ){
				int max_a   = -1;
				float max_s = -100000.0f;
				for( int a=0; a<real_aadim; ++a ){
					float arg = col_freqs[i][a] / freq[a];
					if( arg>max_s ){
						max_s = arg;
						max_a = a;
					}
				}
				//is is an x column?
				if( max_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[max_a];

				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
			
			break;
		}
		case 3:		// min. rel. entropy between p(a) and P(a|x)
		{
			const float** const mat = matrix->get_matrix();
			for( size_t i=0; i<columns; ++i ){
				int max_a   = -1;
				float max_s = -100000.0f;
				for( int a=0; a<real_aadim; ++a ){
					float sum = 0.0f;
					for( int x=0; x<real_aadim; ++x ){
						sum += col_freqs[i][x] * mat[x][a];
					}
					if( sum>max_s ){
						max_s = sum;
						max_a = a;
					}
				}
						
				//is is an x column?
				if( max_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[max_a];

				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
			
			break;
		}
		case 4:		// min. quadratic diff. between profile score and matrix score
		{
			const float * const freq = matrix->get_p_back();
			const float** const mat = matrix->get_matrix();
			for( size_t i=0; i<columns; ++i ){
				int min_a   = -1;
				float min_s = 100000.0f;
				for( int a=0; a<real_aadim; ++a ){
					float sum = 0.0f;
					for( int x=0; x<real_aadim; ++x ){
						sum += col_freqs[i][x] * pow( (fast_log2(col_freqs[i][x] / freq[x]) - (0.5 * mat[x][a])) , 2);
					}
					if( sum<min_s ){
						min_s = sum;
						min_a = a;
					}
				}
				//is is an x column?
				if( min_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[min_a];
				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
			
			break;
		}	
		case 5:		// min. expected score difference between homologs and random seqs
		{
			const float * const freq = matrix->get_p_back();
			const float** const mat = matrix->get_matrix();
			for( size_t i=0; i<columns; ++i ){
				int max_a   = -1;
				float max_s = -100000.0f;
				for( int a=0; a<real_aadim; ++a ){
					float sum = 0.0f;
					for( int x=0; x<real_aadim; ++x ){
						sum += (col_freqs[i][x] - freq[x]) * mat[x][a];
					}
					if( sum>max_s ){
						max_s = sum;
						max_a = a;
					}
				}
				//is is an x column?
				if( max_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[max_a];

				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
			
			break;
		}
		case 6:		// max. profile score minus background frequency
		{
			const float * const freq = matrix->get_p_back();
			for( size_t i=0; i<columns; ++i ){
				int max_a   = -1;
				float max_s = -100000.0f;
				for( int a=0; a<real_aadim; ++a ){
					float arg = col_freqs[i][a] - freq[a];
					if( arg>max_s ){
						max_s = arg;
						max_a = a;
					}
				}
				//is is an x column?
				if( max_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[max_a];

				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
			
			break;
		}			
		case 1:		// conventional method
		default:
		{
			for( size_t i=0; i<columns; ++i ){
				int max_a   = -1;
				float max_s = 0.0f;
				for( int a=0; a<real_aadim; ++a ){
					if( col_freqs[i][a]>max_s ){
						max_s = col_freqs[i][a];
						max_a = a;
					}
				}
				//is is an x column?
				if( max_a==-1 ){
					if( col_freqs[i][AMINOACID_DIM-1]>0.0f ){
		 				str << 'X';
					}else{
						std::cerr << "Strange column found at " << i << " (" << filename << ")" << std::endl;
					}
				}else str << int2aa[max_a];

				if(line_width>0 && (i+1)%line_width==0 && i>0) str << std::endl;	
			}
			if( line_width>0 && columns%line_width!=0 ) str << std::endl;
	
			break;
		}
	}
	
	for( size_t i=0; i<columns; ++i){
		delete [] col_freqs[i];
	}
	delete [] col_freqs;
	return str.str();
}

/////////////////////////////////////////////////////////////////////////////////////
// fast log base 2
/////////////////////////////////////////////////////////////////////////////////////

// This function returns log2 with a max abolute deviation of +/- 1.5E-5 (typically 0.8E-5). 
// It takes 1.42E-8 s  whereas log2(x) takes 9.5E-7 s. It is hence 9.4 times faster. 
// It makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign, 
// the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
// The following 23 bits, m, give the mantisse, the binary digits behind the decimal point.
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127) 
// The expression (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee, i.e. 
// the largest integer that is smaller than log2(x) (e.g. -1 for 0.9). *(int *)&x is an integer which 
// contains the bytes as the floating point variable x is represented in memory.
// Check:  assert( sizeof(f) == sizeof(int) );
// Check:  assert( sizeof(f) == 4 );
inline float Alignment::fast_log2(float x) 
{
  static float lg2[1025];         // lg2[i] = log2[1+x/1024]
  static float diff[1025];        // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
  static char initialized;
  if (x<=0) return -100000;
  if (!initialized)   //First fill in the arrays lg2[i] and diff[i]
    {
      float prev = 0.0f;
      lg2[0] = 0.0f;
      for (int i=1; i<=1024; ++i) 
        {
          lg2[i] = log(float(1024+i))*1.442695041-10.0f;
          diff[i-1] = (lg2[i]-prev)*1.2352E-4;
          prev = lg2[i];
        }
      initialized=1;
    }  
  int a = (((*((int *)&x)) & 0x7F800000) >>23 )-0x7f;
  int b =  ((*((int *)&x)) & 0x007FE000) >>13;
  int c =  ((*((int *)&x)) & 0x00001FFF);
  return a + lg2[b] + diff[b]*(float)(c);
}
