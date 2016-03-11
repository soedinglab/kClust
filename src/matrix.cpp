/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "matrix.h"


Matrix::Matrix(const Matrix::mtype t)throw (std::exception):AMINOACID_DIM(21) {
	matrix = new float*[AMINOACID_DIM];
	for( size_t i=0; i<AMINOACID_DIM; ++i ) matrix[i] = new float[AMINOACID_DIM];
	
	diagonal_scores = new float[AMINOACID_DIM];

	original = new float*[AMINOACID_DIM];
	for( size_t i=0; i<AMINOACID_DIM; ++i ) original[i] = new float[AMINOACID_DIM];
	
	// background aminoacid probabilities
	p_background = new float[AMINOACID_DIM];
	for( size_t i=0; i<AMINOACID_DIM; ++i ) p_background[i] = 0.0f;
	
	p_query = new float[AMINOACID_DIM];

	int2aa  = new char[AMINOACID_DIM];
	int2aa[AMINOACID_DIM-1] = 'X';
	
	// set x-row and -column 
	for( size_t i=0; i<AMINOACID_DIM; ++i){
		matrix[i][AMINOACID_DIM-1] = -1.0f;
		matrix[AMINOACID_DIM-1][i] = -1.0f;
		original[i][AMINOACID_DIM-1] = 0.0f;
		original[AMINOACID_DIM-1][i] = 0.0f;
	}
	
	const size_t char_size = (int)pow(2, 8*sizeof(char));
	aa2int = new int[char_size];
	for(size_t i=0; i<char_size; ++i) aa2int[i]=-1;
	
	if( t!=static_blosum62 ) _read_blosum_matrix(_get_fn(t).c_str());
	else             _copy_dummy();

	for(size_t i=0; i<AMINOACID_DIM; ++i) aa2int[int2aa[i]] = i;
	
	aa_powers = new size_t[32];
	size_t pow = 1;
	for( size_t i=0; i<32; ++i ){
		aa_powers[i] = pow;
		pow *= AMINOACID_DIM;
	}
	// caluculate matrix entries according to the formula matrix[i][j] = 1/scale * log (p(i,j)/(p(i), p(j)))
	reset();
	
	_dimers = new dimer*[aa_powers[2]];
	for( size_t i=0; i<aa_powers[2]; ++i ){
		_dimers[i] = new dimer[aa_powers[2]];	
	}
	_init_dimers();
	_init_main_diagonal_scores();

	//_print_debug(std::cerr);
	negatives=0;
	error_sum=0.0;
	
}


Matrix::~Matrix(){
	for( size_t i=0; i<AMINOACID_DIM; ++i ) delete [] matrix[i];
	delete [] matrix;
	for( size_t i=0; i<AMINOACID_DIM; ++i ) delete [] original[i];
	delete [] original;
	delete [] p_background;
	delete [] diagonal_scores;
	for (size_t i = 0; i < aa_powers[2]; i++) delete [] _dimers[i];
	delete [] _dimers;
	delete [] p_query;
	delete [] int2aa;
	delete [] aa2int;
	delete [] aa_powers;
}

void Matrix::_init_main_diagonal_scores(){
//	std::cout << "max scores: ";
	for( size_t a=0; a<AMINOACID_DIM; ++a ) {
		diagonal_scores[a] = matrix[a][a];
//		std::cout << matrix[a][a] << " ";
	}
//	std::cout << "\n";
}

void Matrix::_copy_dummy(){
	for(size_t i=0; i<21; ++i)
		for(size_t j=0; j<21; ++j)
			original[i][j] = __DUMMY_MATRIX[i][j];
	for(size_t i=0; i<21; ++i) int2aa[i]       = __DUMMY_INT2AA[i];
	for(size_t i=0; i<21; ++i) p_background[i] = __DUMMY_P_BACKGROUND[i];
	scale = __DUMMY_SCALE;
}

void Matrix::init_p_query( Sequence *s ){
	const char *seq_ptr = s->get_sequence();
	const double L      = s->length();
	const double pseudo_counts = 100.0;
	for( size_t i=0; i<20; ++i) p_query[i] = pseudo_counts * p_background[i];
	for( size_t i=0; i<L; ++i) ++p_query[seq_ptr[i]];	
	double sum=0.0;
	for( size_t i=0; i<20; ++i) sum += p_query[i];
	for( size_t i=0; i<20; ++i) p_query[i]/=sum;
}

const float** const Matrix::get_matrix()           const { return (const float**)matrix;   }
const float* const Matrix::get_diagonal_scores()   const { return (const float*) diagonal_scores; }
const float** const Matrix::get_prob_matrix()      const { return (const float**)original; }
const float* const Matrix::get_p_query()           const { return (const float*)p_query;   }
const float* const Matrix::get_p_back()            const { return (const float*)p_background;   }
const size_t Matrix::get_aa_dim()                  const { return AMINOACID_DIM;           }
const char* const Matrix::get_int2aa()             const { return int2aa;                  }
const int* const Matrix::get_aa2int()              const { return aa2int;                  }


// Rescale BLOSUM matrix for query amino acid frequencies with conjugate gradient method. Used in kDP
void Matrix::wn_rescale(Sequence *s) throw (std::exception){
	init_p_query(s);
	parameter_wrapper pwrap;
	pwrap.matrix  = original;
	pwrap.p_query = p_query;
	int code=0;
  	double val_min;
  	const int dim=20;
  	const int max_iterations=40;
	double x[dim];
	for( size_t a=0; a<20; ++a) x[a] = pow(p_query[a]/p_background[a], 0.7121);
	wn_conj_gradient_method(&code, &val_min, x, dim, &_func, &_grad, max_iterations, &pwrap);
	if( code!=WN_SUCCESS && code!=WN_SUBOPTIMAL) 
		throw MyException("wn_conj_gradient_method failed: return-code is '%i'\n", code);
	for (int a=0; a<20; ++a)
		for (int b=0; b<20; ++b) 
		matrix[a][b] = 2.0*_log2(x[a]*original[a][b]*x[b]/p_query[a]/p_query[b]);
	bool bprint=false;
	double za[20] = {0.0f};
	for (int a=0; a<20; ++a){
		if( x[a]<0.0 ) bprint=true;
		for (int b=0; b<20; ++b) za[a] += original[a][b] * x[b];
	}
	double error=0.0;
	for (int a=0; a<20; ++a){
		double tmp = 1.0-(za[a]*x[a])/p_query[a];
		error += tmp*tmp;
	}
	error_sum += error;
	if(bprint) 
		throw MyException("One or more negative scaling factor(s) occurred while rescaling the amino acid substitution matrix for '%s'!", s->get_header());
	
	_init_main_diagonal_scores();
/*
	if(bprint || error >0.01){
		if(bprint) ++negatives;
		printf("Error: %5.5f\n", error);
		fprintf(stderr, "Minimizing factors (x[])\n");
		for( int i=0; i<20; ++i){
			fprintf(stderr, "%2.4f ", x[i]);
		}
		fprintf(stderr, "\n");
		for (int a=0; a<20; ++a){
			fprintf(stderr, "%2.4f ", (za[a]*x[a])/p_query[a]);
		}
		fprintf(stderr, "\n");
		
	}*/
}

double _func(double *x, void *params){	
	Matrix::parameter_wrapper *p = (Matrix::parameter_wrapper*)params;
	const float** const matrix   = (const float** const)p->matrix;
	const float* const p_query   = p->p_query;
	double sum=0.0;
	for (int a=0; a<20; ++a){
		double za = 0.0;
		for (int b=0; b<20; ++b) za += matrix[a][b] * x[b];
		za   = (za*x[a]) - p_query[a];
		sum += za*za;
	}
	return sum;
}

//the gradient (grad[20]) of the function at the point x[20]
void _grad(double *grad, double *x, void *params){
	Matrix::parameter_wrapper *p = (Matrix::parameter_wrapper*)params;
	const float** const matrix = (const float** const)p->matrix;
	const float* const p_query = p->p_query;
	double za[20] = {0.0f};
	for (int a=0; a<20; ++a)
		for (int b=0; b<20; ++b) za[a] += matrix[a][b] * x[b];
	double wa[20];
	for (int a=0; a<20; ++a) wa[a] = (x[a]*za[a]) - p_query[a];
	double tmp;
	for (int a=0; a<20; ++a){
		tmp=0.0;
		for (int b=0; b<20; ++b) tmp += wa[b]*x[b]*matrix[b][a];
		tmp += wa[a]*za[a];
		tmp *= 2.0;
		grad[a] = tmp;
	}
	return;
}

void Matrix::reset(){
	for(size_t i=0; i<20; ++i){
		for(size_t j=0; j<20; ++j){
			matrix[i][j] = _log2(original[i][j]/(p_background[i]*p_background[j])) / scale;
//			std::cout << matrix[i][j] << "  ";
		}
//		std::cout << matrix[i][20] << "\n";
	}
/*	std::cout << "\n";
	for(size_t j=0; j<21; ++j){
		std::cout << matrix[20][j] << "  ";
	}
	std::cout << "\n";*/
	_init_main_diagonal_scores();
}

void Matrix::round_bit_scores(){
	for(size_t i=0; i<20; ++i)
		for(size_t j=0; j<20; ++j)
			matrix[i][j] = round( matrix[i][j]);
	_init_main_diagonal_scores();
}

void Matrix::_copy(){
	for(size_t i=0; i<21; ++i)
		for(size_t j=0; j<21; ++j)
			matrix[i][j] = original[i][j];
	_init_main_diagonal_scores();
}

std::ostream& Matrix::print(std::ostream &out){
	out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
	char buf[20];
	out << " ";
	for(size_t i=0; i<21; ++i){
		sprintf(buf, "   %c  ", int2aa[i]);
		out << buf;
	}
	out << std::endl;
	for(size_t i=0; i<21; ++i){
		for(size_t j=0; j<21; ++j){
			sprintf(buf, "%+2.2f ", matrix[i][j]);
			out << buf;
		}
		out << std::endl;
	}
	return out;
}

std::ostream& Matrix::print_frequencies(std::ostream &out){
	out << "{" << std::endl;
	for(size_t i=0; i<21; ++i){
		out << "{ ";
		for(size_t j=0; j<21; ++j){
			out << original[i][j];
			if( j!=20 ) out << ", ";
			else out << "}";
		}
		if(i!=20) out << "," << std::endl;
		else      out << "}" << std::endl;
	}
	return out;
}


std::ostream& Matrix::print_p_query(std::ostream &out){
	for(size_t j=0; j<21; ++j)
		out << int2aa[j] << " : " << p_query[j] << " " <<  p_background[j] << std::endl;
	return out;
}

std::ostream& Matrix::print_p_back(std::ostream &out){
	for(size_t j=0; j<21; ++j)
		out << j << " " << int2aa[j] << " " << p_background[j] << std::endl;
	return out;
}


std::ostream& Matrix::print_int(std::ostream &out){
	out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
	char buf[20];
	out << " ";
	for(size_t i=0; i<21; ++i){
		sprintf(buf, "   %c", int2aa[i]);
		out << buf;
	}
	out << std::endl;
	for(size_t i=0; i<21; ++i){
		sprintf(buf, "%c ", int2aa[i]);
		out << buf;
		for(size_t j=0; j<21; ++j){
			sprintf(buf, "%+3i ", (int)round(matrix[i][j]));
			out << buf;
		}
		out << std::endl;
	}
	return out;
}

const std::string Matrix::_get_fn(const Matrix::mtype t){
	std::string p="./";
	if(getenv("BLOSUM_MATRICES")) p = std::string(getenv("BLOSUM_MATRICES"));
	switch(t){
		case Matrix::blosum30: return (p+"/blosum30.out");
		case Matrix::blosum35: return (p+"/blosum35.out");
		case Matrix::blosum40: return (p+"/blosum40.out");
		case Matrix::blosum45: return (p+"/blosum45.out");
		case Matrix::blosum50: return (p+"/blosum50.out");
		case Matrix::blosum55: return (p+"/blosum55.out");
		case Matrix::blosum60: return (p+"/blosum60.out");
		case Matrix::blosum62: return (p+"/blosum62.out");
		case Matrix::blosum65: return (p+"/blosum65.out");
		case Matrix::blosum70: return (p+"/blosum70.out");
		case Matrix::blosum75: return (p+"/blosum75.out");
		case Matrix::blosum80: return (p+"/blosum80.out");
		case Matrix::blosum85: return (p+"/blosum85.out");
		case Matrix::blosum90: return (p+"/blosum90.out");
		case Matrix::blosum95: return (p+"/blosum95.out");
		case Matrix::blosum100: return (p+"/blosum100.out");
		case Matrix::static_blosum62 : return ("no-file");
	}
	return "";
}

void Matrix::_read_blosum_matrix(const char *fn) throw (std::exception){
	std::ifstream in(fn);
	if( in.fail() ) 
		throw MyException("Cannot read '%s'!\nPlease export $BLOSUM_MATRICES=/dir-where-the-blosumXX.out-files-reside/", fn);
	int c      = 0;
	int row    = 0;
	int column = 0;
	std::string line;
	bool capture = false;
	while( in.good() ){
		getline( in, line );
		if( line.length()>11 && line.substr(0, 11)!="Frequencies" && !capture ) continue;
		if( line.length()>11 && line.substr(0, 11)=="Frequencies"){
			capture=true;
			continue;
		}
		if( row==20 ) break;
		std::stringstream stream(line); std::string h; stream >> h;
		if( h=="" ) continue;
		if( isalpha(h.at(0)) ){
			int2aa[c++] = toupper( h.at(0) );
			while(	stream >> h ){
				int2aa[c++] = toupper( h.at(0) );
				if( c>20 ) throw MyException("Blosum matrix file '%s' has wrong format!\n", fn);
			}
		}else{
			column = 0;
			stream.clear();
			stream.str(line);
			float f;
			while( stream >> f ){
				original[row][column] = f;
				original[column][row] = f;
				++column;
			}
			++row;
		}
	}
	if( c!=20 ) throw MyException("Blosum matrix file '%s' has wrong format!\n", fn);
	in.close();

	float sum=0.0f;
	for(size_t i=0; i<20; ++i)
		for(size_t j=0; j<20; ++j){
			if( i==j ) p_background[i] += original[i][j];
			else       p_background[i] += (original[i][j]/2.0f);
			if( j<=i ) sum += original[i][j]; 
		}

	const float _2sum = 2.0*sum;	
	float pbsum = 0.0f;
	for(size_t i=0; i<20; ++i){
		pbsum += p_background[i];
		for(size_t j=0; j<20; ++j)
			if( i==j ) original[i][j] = original[i][j] / sum;
			else       original[i][j] = original[i][j] / _2sum;
	}

	for(size_t i=0; i<20; ++i)p_background[i] /= sum;

	float entropy=0.0f;
	//compute entropy	
	for(size_t i=0; i<20; ++i)
		for(size_t j=0; j<=i; ++j)
			entropy += original[i][j] * _log2( original[i][j] / (p_background[i]*p_background[j]) ) ;
	//set scaling factor for blosum matrices half-bits, third-bits,...
	//std::cerr << "Entropy:" << entropy << std::endl;
	scale = std::min( 0.5, 1.0/(round(2.0/sqrt(entropy))) );
}

std::ostream& Matrix::_print_debug(std::ostream &out){
	char buf[20];
	out << "Frequencies of amino acid substitutions (unscaled, float precision):" << std::endl;
	for(size_t i=0; i<21; ++i){
		sprintf(buf, "%6c ", int2aa[i]);
		out << buf;
	}
	out << std::endl;
	for(size_t i=0; i<21; ++i){
		for(size_t j=0; j<21; ++j){
			sprintf(buf, "%+2.3f ", original[i][j]);
			out << buf;
		}
		out << std::endl;
	}
	out << std::endl;
	out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
	for(size_t i=0; i<21; ++i){
		sprintf(buf, "%6c ", int2aa[i]);
		out << buf;
	}
	out << std::endl;
	for(size_t i=0; i<21; ++i){
		for(size_t j=0; j<21; ++j){
			sprintf(buf, "%+2.3f ", matrix[i][j]);
			out << buf;
		}
		out << std::endl;
	}
	out << std::endl;
	return out;
}	


void Matrix::read_matrix(const float src[][21]){
	for(size_t i=0; i<21; ++i)
		for(size_t j=0; j<21; ++j)
			matrix[i][j] = src[i][j];
	_init_main_diagonal_scores();
}

void Matrix::_init_dimers(){
	const size_t aa_2 = AMINOACID_DIM*AMINOACID_DIM;
	dimer **to_sort = new dimer*[aa_2];
	for( size_t i=0; i<aa_2; ++i ) to_sort[i] = new dimer();
	for( size_t i=0; i<AMINOACID_DIM; ++i ){
		for( size_t j=0; j<AMINOACID_DIM; ++j ){
			// unique index of the dimer with base 2
			size_t at = i + j*AMINOACID_DIM;
			for( size_t p1=0; p1<AMINOACID_DIM; ++p1 ){
				for( size_t p2=0; p2<AMINOACID_DIM; ++p2 ){
					float score = matrix[i][p1] + matrix[j][p2];
					size_t idx  = p1+ p2*AMINOACID_DIM ;
					to_sort[idx]->score = score;
					to_sort[idx]->idx   = idx;
				}
			}		
			dimer **sorted = _merge_sort_dimers( to_sort, aa_2 );
			for( size_t t=0; t<aa_2; ++t ){
				_dimers[at][t].idx   = sorted[t]->idx;
				_dimers[at][t].score = sorted[t]->score;
			}
			delete [] sorted;
		}	
	}
	/*
	//print
	for( size_t i=0; i<AMINOACID_DIM; ++i ){
		for( size_t j=0; j<AMINOACID_DIM; ++j ){
			size_t at = i*AMINOACID_DIM + j;
			std::cout <<int2aa[i]<<int2aa[j]<<":";
			for( size_t p=0; p<aa_powers[2]; ++p ){
				size_t idx  = _dimers[at][p].idx;
				size_t a2 = idx%AMINOACID_DIM;
				size_t a1 = (idx-a2)/AMINOACID_DIM;
				float score = _dimers[at][p].score;
				std::cout<<int2aa[a1]<<int2aa[a2]<< " " << _dimers[at][p].idx << " ("<<score<<") ";
			}
			std::cout << std::endl;
		}	
	}
	*/
	for( size_t i=0; i<aa_2; ++i ) delete to_sort[i];
	delete [] to_sort;
}

Matrix::dimer** Matrix::_merge_sort_dimers( dimer **dimers, size_t len ){
	if( len<2 ) return dimers;
	size_t l=0;
	size_t h=0;
	size_t m=(len-1)/2;	
	Matrix::dimer **lower = _merge_sort_dimers( dimers, m+1 );
	Matrix::dimer **upper = _merge_sort_dimers( dimers+m+1, len-m-1 );
	Matrix::dimer **tmp = new dimer*[len];
	while (l+h<len){ 
		if( l<=m && (h>=len-m-1 || lower[l]->score > upper[h]->score ) ){
			tmp[l+h] = lower[l];	
			++l;
		}else{ 
			tmp[l+h] = upper[h];
			++h;
		}
	}
  	if (len>2) delete [] lower; 
  	if (len>3) delete [] upper; 
  	return tmp;
}
