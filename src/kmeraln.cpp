/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "kmeraln.h"

Kmeraln::Kmeraln(	const int k, 
			Matrix *mat, 
			const float p_m,
			const int delta,
			const float G, 
			const float E, 
			const float F,
			const float Lamda, 
			const float K, 
			const float H,
			const double db_len,
			const float sum_log_db_len, 
			const float db_seqs, 
			const int KMER_OFFSET,
			const float COV ) throw (std::exception):
			AMINOACID_DIM( mat->get_aa_dim() ),
			int2aa( mat->get_int2aa() ),
			p_m(p_m),
			delta(delta),
			G(G), 
			E(E), 
			F(F),
			Lamda(Lamda),
			K(K),
			H(H),
			db_len(db_len),
			sum_log_db_len(sum_log_db_len), 
			db_seqs(db_seqs), 
			_dropoff(40.0f),
			_grid_score(0.0f),
			_k(k), 
			COV(COV) {
	ExtMatchHeap  = new HeapWrapper<_extended_match>( ((4*1024*1024)-32)/sizeof(_extended_match) );
	_trace        = 0;
	_score        = 0.0f;
	aln_cov_short = 0.0f;
	aln_cov_long  = 0.0f;
	_first        = 0;
	forward       = new float[32768];
	backward      = new float[32768];
	if( KMER_OFFSET<0 )  _kmer_offset =  (_k-1)/2;
	else                 _kmer_offset    = KMER_OFFSET;
	aln = 0;
	
}

Kmeraln::~Kmeraln(){
	delete ExtMatchHeap;
	delete [] forward;
	delete [] backward;
	delete aln;
}

//x query sequence (the sequence for which the similar k-mers were created)
//y the representantive sequence in the table
void Kmeraln::align(	Kmeraln::Method method, 
			Recycle_table::_match *mlist, 
			Sequence *x, 
			Sequence *y) throw (std::exception) {
	this->Lx      = x->length();
	this->Ly      = y->length();
	this->X_seq   = x;
	this->Y_seq   = y;
	this->X       = x->sequence();
	this->Y       = y->sequence();
	_trace        = 0;
	_score        = 0.0f;
	aln_cov_short = 0.0f;
	aln_cov_long  = 0.0f;
	_first        = 0;
//	std::cerr << "Sequence X (query): " << x->get_header() << "\n";
//	std::cerr << "Sequence Y (representative): " << y->get_header() << "\n";
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "   Length of sequence X                 : " << Lx << std::endl;
		std::cerr << "   Length of sequence Y                 : " << Ly << std::endl;
	#endif

	delete aln;
	aln=0;
	_init_list(mlist);
	switch(method){
		case Kmeraln::FAST_ADDR:{
			_init_fast_addr();
			_fill_fast();
			_delete_fast_addr();
			aln = traceback();
			break;}
		case Kmeraln::FULL_BAND:{
			_fill();
			aln = traceback();
			break;}
		default: throw MyException("Unknown method for kmeraln!");
	}
}

size_t Kmeraln::get_memory_usage(){
	size_t mem=0;
	mem += 32768*2*sizeof(float);
	mem += ExtMatchHeap->get_instance_memory_usage();
	return mem;
}

std::ostream &Kmeraln::print(	std::ostream& out,
										Sequence *x, 
										Sequence *y, 
										const float** const matrix, 
										const size_t width) throw (std::exception){
	if( aln==0 ) throw MyException("No alignment available!");
	out << std::endl;
	out << "Sequence x(upper)                  : " << x->get_header()           << std::endl;
	out << "Sequence y(lower)                  : " << y->get_header()           << std::endl;
	out << "Length of sequence x               : " << aln->len_x            << std::endl;
	out << "Number of aligned residues of x    : " << aln->aln_len_x        << std::endl;
	out << "Length of sequence y               : " << aln->len_y            << std::endl;
	out << "Number of aligned residues of y    : " << aln->aln_len_y        << std::endl;
	out << "Alignment start position in x      : " << aln->x_aln_start      << std::endl;
	out << "Alignment end position in x        : " << aln->x_aln_end        << std::endl;
	out << "Alignment start position in y      : " << aln->y_aln_start      << std::endl;
	out << "Alignment end position in y        : " << aln->y_aln_end        << std::endl;
	out << "Length of alignment                : " << aln->len_aln          << std::endl;
	out << "Score                              : " << aln->score            << std::endl;
	out << "Identities                         : " << aln->idents           << std::endl;
	out << "Positive matches (incl. idents)    : " << aln->pos_matches      << std::endl;
	out << "Less or equal to zero matches      : " << aln->neg_zero_matches << std::endl;
	out << "Gaps                               : " << aln->gaps             << std::endl;	
	out << std::endl;
	out << "Normalized score:                  : " << _score/std::min(Lx,Ly)<< std::endl;
	out << "Sequence identity                  : " << aln->idents/(float)std::min(Lx,Ly)  << std::endl;
	out << "Alignment coverage short           : " << aln_cov_short                       << std::endl;
	out << "Alignment coverage long            : " << aln_cov_long                        << std::endl;
	out << std::endl;

	char *aln_seq_x = new char[Lx+Ly];
	char *aln_seq_y = new char[Lx+Ly];

	const size_t len_aln = aln->len_aln;

	for( size_t i=0; i<len_aln; ++i){
		if( aln->x[i]==-1 ){
			aln_seq_x[i] = '-';
			aln_seq_y[i] = Y[ aln->y[i] ];
		}else if( aln->y[i]==-1 ){
			aln_seq_x[i] = X[ aln->x[i] ];
			aln_seq_y[i] = '-';
		}else{
			aln_seq_x[i] = X[ aln->x[i] ];
			aln_seq_y[i] = Y[ aln->y[i] ];
		}
	}
	
	char buffer[width];
	char xbuffer[width];
	char ybuffer[width];

	size_t jj, i, j;

	for(i=0; i<len_aln; i+=width ){
		for(jj=0, j=i; j<std::min(i+width, len_aln); ++j){
			if( aln_seq_x[j]==aln_seq_y[j] && aln_seq_x[j]!='-' )
				buffer[jj] = '|';
			else if( (aln_seq_x[j]!='-') && (aln_seq_y[j]!='-')){
				if( (matrix[aln_seq_x[j]][aln_seq_y[j]] > 0) )
					buffer[jj] = ':';
				else
					buffer[jj] = '~';
			}else buffer[jj] = ' ';
			if( aln_seq_x[j]<AMINOACID_DIM ) xbuffer[jj] = int2aa[ aln_seq_x[j] ];
			else xbuffer[jj] = aln_seq_x[j];
			if( aln_seq_y[j]<AMINOACID_DIM ) ybuffer[jj] = int2aa[ aln_seq_y[j] ];
			else ybuffer[jj] = aln_seq_y[j];
			++jj;
		}
		out.write(xbuffer, jj);	
		out << std::endl;
		out.write(buffer, jj);
		out << std::endl;
		out.write(ybuffer, jj);
		out << std::endl << std::endl;
	}
	out << "--------------------------------------------------------"           << std::endl;
//	#ifdef __CM_DEBUG_KMERALN
//		_compute_overlap_score(aln_seq_x, aln_seq_y, aln->len_aln, matrix);
//	#endif
	delete [] aln_seq_x;
	delete [] aln_seq_y;
	return out;
}

std::ostream &Kmeraln::print_real_trace(std::ostream& out) throw (std::exception){

	if( aln==0 ) throw MyException("No alignment available!");
	
	char *seq_aln_x = new char[Lx+Ly];
	char *seq_aln_y = new char[Lx+Ly];
	const size_t len_aln = aln->len_aln;

	for( size_t i=0; i<len_aln; ++i){
		if( aln->x[i]==-1 ){
			seq_aln_x[i] = '-';
			seq_aln_y[i] = Y[ aln->y[i] ];
		}else if( aln->y[i]==-1 ){
			seq_aln_x[i] = X[ aln->x[i] ];
			seq_aln_y[i] = '-';
		}else{
			seq_aln_x[i] = X[ aln->x[i] ];
			seq_aln_y[i] = Y[ aln->y[i] ];
		}
	}
	
	size_t l=aln->x[0], k=aln->y[0];
	for( size_t i=0; i<len_aln; ++i ){
		out << l << " " << k << std::endl;
		if( seq_aln_x[i]!='-' ) ++l;
		if( seq_aln_y[i]!='-' ) ++k;
	}
	
	delete [] seq_aln_x;
	delete [] seq_aln_y;

	return out;
}

std::ostream &Kmeraln::print_alignment(std::ostream& out) throw (std::exception){

	if( aln==0 ) throw MyException("No alignment available!");
	
	char *seq_aln_x = new char[Lx+Ly];
	char *seq_aln_y = new char[Lx+Ly];
	const size_t len_aln = aln->len_aln;

	int matches = 0;
	int start_pos_x = -1;
	int end_pos_x = -1;
	int start_pos_y = -1;
	int end_pos_y = -1;
	
	for( size_t i=0; i<len_aln; ++i){
		if( aln->x[i]==-1 ){
			seq_aln_x[i] = '-';
			seq_aln_y[i] = int2aa[Y[ aln->y[i] ]];
			if (start_pos_y == -1) start_pos_y = aln->y[i];
			end_pos_y = aln->y[i];
		}else if( aln->y[i]==-1 ){
			seq_aln_x[i] = int2aa[X[ aln->x[i] ]];
			seq_aln_y[i] = '-';
			if (start_pos_x == -1) start_pos_x = aln->x[i];
			end_pos_x = aln->x[i];
		}else{
			seq_aln_x[i] = int2aa[X[ aln->x[i] ]];
			seq_aln_y[i] = int2aa[Y[ aln->y[i] ]];
			if (start_pos_x == -1) start_pos_x = aln->x[i];
			if (start_pos_y == -1) start_pos_y = aln->y[i];
			end_pos_x = aln->x[i];
			end_pos_y = aln->y[i];
			if (X[ aln->x[i] ] == Y[ aln->y[i] ]) matches++;
		}
	}
	
	float seq_id = (float)matches/(float)len_aln;
	
	size_t l=aln->x[0], k=aln->y[0];
	
	out << "#" << seq_id << "|" << _score << "|" << (start_pos_x+1) << ":" << (end_pos_x+1) << "|" << (start_pos_y+1) << ":" << (end_pos_y+1) << "\n";
	out << X_seq->get_header() << "\n";
	for(int i = 0; i < len_aln; ++i ){
		out << seq_aln_x[i];
	}
	out << "\n";
	out << Y_seq->get_header() << "\n";
	for(int j = 0; j < len_aln; ++j ){
		out << seq_aln_y[j];
	}
	out << "\n";
	
	delete [] seq_aln_x;
	delete [] seq_aln_y;

	return out;
}

void Kmeraln::_init_list(Recycle_table::_match *match) throw (std::exception){

	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "Initialization of match-list ..." << std::endl;
	#endif

	if( match==0 ){
		_first=0;
		return; 
	}
	
	int d_min=0, d_max=0, d, skipped;
	const int Lmax = std::max(Lx, Ly);
	const int Lmin = std::min(Lx, Ly);
	if( COV>0.0 && (Lmin/(double)Lmax) <= COV ) 
		throw MyException("The coverage criterion (%2.2f) is not fulfilled with these sequences (x:%i y:%i)!", COV, Lx, Ly); 
	if(COV>0.0f){
		d_min   = (int)(floor(COV*Lmax) - Ly);
		d_max   = (int)(Lx - floor(COV*Lmax));
		#ifdef __CM_DEBUG_KMERALN
		std::cerr << "   Restriction to band -> d_min:" << d_min << " d_max:" << d_max << std::endl;
		#endif
		skipped = 0;
	}
	
	ExtMatchHeap->reset();
	_extended_match *tmp      = ExtMatchHeap->get_next();
	_extended_match *previous = tmp;
	tmp->i         = 0;
	tmp->j         = 0;
	tmp->kmerscore = _grid_score;
	tmp->dummy     = true;
	_first         = tmp;

	// init synthetic matches
	int trigger_synthetic_matches_i = 0;
	for( int j=delta; j<Ly; j+=delta ){
		#ifdef __CM_DEBUG_KMERALN_MATCH_LIST
			std::cerr << "   Adding artificial match with properties: i=";
			std::cerr << trigger_synthetic_matches_i;
			std::cerr << "      j=" << j << " score=" << 0.0f << std::endl;
		#endif
		if( COV>0.0f ){
			if( -j<=d_min || -j>=d_max ) continue;
		}
		tmp              = ExtMatchHeap->get_next();
		tmp->i           = 0;
		tmp->j           = j;
		tmp->kmerscore   = _grid_score;
		tmp->dummy       = true;
		previous->inext  = tmp; 
		previous         = tmp;
	}
	trigger_synthetic_matches_i += delta;

	int matches_num = 0;
	while( match!=0 ){
		while( match->i >= trigger_synthetic_matches_i ){
			for( int j=0; j<Ly; j+=delta ){
//				if ((match->i == trigger_synthetic_matches_i) && (match->j == j)) continue;
				#ifdef __CM_DEBUG_KMERALN_MATCH_LIST
					std::cerr << "   Adding artificial match: (" << trigger_synthetic_matches_i << "," << j << ")\n";
				#endif
				if( COV>0.0f ){
					d               = trigger_synthetic_matches_i-j;
					if( d<=d_min || d>=d_max ) continue;
				}
				tmp             = ExtMatchHeap->get_next();
				tmp->i          = trigger_synthetic_matches_i;
				tmp->j          = j;
				tmp->kmerscore  = _grid_score;
				tmp->dummy      = true;
				previous->inext = tmp; 
				previous        = tmp;
			}
			trigger_synthetic_matches_i += delta;
		}
		#ifdef __CM_DEBUG_KMERALN_MATCH_LIST
			std::cerr << "   Copying match with properties: (" << match->i << "," << match->j << ") score=" << match->score << std::endl;
		#endif
		if( COV>0.0f ){
			d = match->i - match->j;
			if( d<=d_min || d>=d_max ){ 
				match = match->next;
				++skipped;
				++matches_num;
				continue;
			}
		}
		tmp              = ExtMatchHeap->get_next();
		tmp->i           = match->i;
		tmp->j           = match->j;
		tmp->kmerscore   = match->score;
		tmp->dummy       = false;
		previous->inext  = tmp; 
		previous         = tmp;
		match            = match->next;
		++matches_num;
	}
	tmp->inext = 0;
	
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "   Copied ";
		std::cerr << ExtMatchHeap->get_number_of_current_objects_in_use();
		std::cerr << " kmer matches."<< std::endl;
		std::cerr << "Initialization of match-list done." << std::endl; 
	#endif
}

void Kmeraln::_init_fast_addr() throw (std::exception){
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "Initialization of fast addressing array ..." << std::endl;
	#endif
	//number of kmers matching randomly in a band of width delta and length=#diagonal=Li+Lj-1
	C        = std::max( (int)( ceil(delta * p_m * (Ly+Lx-1)) ) , 1);
	W        = std::max( (int)ceil((Ly+Lx-1)/(double)C) , 1);
	diag_ptr = new Kmeraln::_extended_match*[C];

	Kmeraln::_extended_match **m=diag_ptr;
	while( m!=diag_ptr+C )*m++=0;

	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "   Probability for kmer-match by chance : " << p_m << std::endl;
		std::cerr << "   Number of diagonals                  : " << (Ly+Lx-1) << std::endl;
		std::cerr << "   Size of address array (C)            : " << C << std::endl;
		std::cerr << "   Width of diagonal segments (W)       : " << W << std::endl;
		std::cerr << "Initialization of fast addressing array done."  << std::endl;
	#endif
}

void Kmeraln::_delete_fast_addr(){
	delete [] diag_ptr;
}

void Kmeraln::_fill_fast(){
	
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "Dynamic programming on kmer-match-list ..." << std::endl;
	#endif

	#ifdef __CM_DEBUG_KMERALN_DYNAMIC
		size_t visited_basepoints  = 0;
		size_t visited_interpoints = 0;
	#endif

	//offset for diagonals to shift into positive range from 0..max_diag
	const int offset = Ly-1;

	Kmeraln::_extended_match *current          = _first;
	Kmeraln::_extended_match *itail            = _first;
	Kmeraln::_extended_match *ihead            = 0;
	Kmeraln::_extended_match *basepoint_runner = 0;
	Kmeraln::_extended_match *intra_runner     = 0;
	Kmeraln::_extended_match *tmp              = 0;

	//inter-diagonal gap length
	int g;
		
	//intra-diagonal gap length
	int z;

	const float k = _k;

	register int d,dx=0,segment,segmentx;
	float max_score, score;
	//hold score of best cell	
	float overall_max_score = -INFINITY;

	while( current!=0 ){

		d        = current->i - current->j;   //diagonal of 'current' match
		segment  = (d+offset)/W;	          //segment of 'current' match
		segmentx = (d-delta+offset)/W;        //segment of the start basepoint
		if( segmentx<0 ) segmentx=0;
		
		#ifdef __CM_DEBUG_KMERALN_DYNAMIC
			std::cerr << "Iteration for " << current << " ";
			std::cerr << " i:"<< current->i << " j:" << current->j;
			std::cerr << " diagonal+offset:" << d+offset <<  " segment:" << segment;
			std::cerr << " score:" << current->kmerscore << std::endl;
		#endif
		
		basepoint_runner = diag_ptr[segmentx];
		max_score        = -INFINITY;
		// basepoint_runner climbs the current basepoint line
		// intra_runner climbs the current diagonal
		while( basepoint_runner!=0 ){
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				++visited_basepoints;
				std::cerr << "     Checking basepoint i:" << basepoint_runner->i ;
				std::cerr << " j:" << basepoint_runner->j;
				std::cerr << " diagonal+offset:" << basepoint_runner->i-basepoint_runner->j + offset;
				std::cerr << " dynscore:" << basepoint_runner->dynscore;
				std::cerr << " addr:" << basepoint_runner<< std::endl;
				std::cerr << "          inter_d_prev:" << basepoint_runner->inter_d_prev << std::endl;
				std::cerr << "          inter_d_next:" << basepoint_runner->inter_d_next << std::endl;
			#endif
//			std::cout << "basepoint: (" << basepoint_runner->i << "," << basepoint_runner->j << ") "<< basepoint_runner<< "\n";
			// diagonal of the basepoint
			dx = basepoint_runner->i - basepoint_runner->j;
			// basepoint is outside the delta window
			if( basepoint_runner->j > current->j ){
				//std::cerr << "     Checking matches on diagonal" << std::endl;
				intra_runner = basepoint_runner->intra_d_next;
				while( intra_runner!=0 ){
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						++visited_interpoints;
						std::cerr << "     Intra diagonal point:"<< intra_runner << " ";
						std::cerr << "i:" << intra_runner->i << " ";
						std::cerr << "j:" << intra_runner->j << " ";
						std::cerr << intra_runner->dynscore << " d: " <<dx << std::endl;
					#endif
//					std::cout << "intra runner: (" << intra_runner->i << "," << intra_runner->j << ")\n";
					if( _in_window(current, intra_runner) ){
//						std::cout << " TRUE\n";
/*						int ret = false;
					    if (current->i == intra_runner->i && current->j == intra_runner->j) 
					    	ret = true;
						else if (current->dummy && intra_runner->dummy)
							ret = (current->i >= intra_runner->i) && (current->j >= intra_runner->j) && (current->i - delta <= intra_runner->i) && (current->j - delta <= intra_runner->j);
						else
							ret = (current->i > intra_runner->i) && (current->j > intra_runner->j) && (current->i - delta <= intra_runner->i) && (current->j - delta <= intra_runner->j);
						if (ret == false) 
							std::cout << "WARNING: from (" << current->i << "," << current->j << "): " << current->dummy << "  to  ""(" << intra_runner->i << "," << intra_runner->j << "):" << intra_runner->dummy <<"  = TRUE\n";
*/
						dx = intra_runner->i-intra_runner->j;
						g  = abs(d-dx)-1;
						//since basepoint.j<j and i>=basepoint.i => dx<d
						z  = current->j-intra_runner->j;
						//the kmer-matches lie on different diagonals
						score = intra_runner->dynscore - (G+(E*g)) - (F*z);
						#ifdef __CM_DEBUG_KMERALN_DYNAMIC
							std::cerr << "          ->Intra-point is valid! score:" << score << std::endl;
						#endif
						if( score>max_score ){
							#ifdef __CM_DEBUG_KMERALN_DYNAMIC
								std::cerr << "          -->Maximum scoring point!" << std::endl;
							#endif
							max_score = score;
							tmp       = intra_runner;
						}
						break;
					}
					intra_runner = intra_runner->intra_d_next;
				}
			}
			//basepoint comes inside the delta window
			else{
				if( _in_window(current, basepoint_runner) ){
//					std::cout << " TRUE\n";
/*					int ret = false;
				    if (current->i == basepoint_runner->i && current->j == basepoint_runner->j) 
				    	ret = true;
					else if (current->dummy && basepoint_runner->dummy)
						ret = (current->i >= basepoint_runner->i) && (current->j >= basepoint_runner->j) && (current->i - delta <= basepoint_runner->i) && (current->j - delta <= basepoint_runner->j);
					else
						ret = (current->i > basepoint_runner->i) && (current->j > basepoint_runner->j) && (current->i - delta <= basepoint_runner->i) && (current->j - delta <= basepoint_runner->j);
					if (ret == false) 
						std::cout << "WARNING: from (" << current->i << "," << current->j << "): " << current->dummy << "  to  ""(" << basepoint_runner->i << "," << basepoint_runner->j << "):" << basepoint_runner->dummy <<"  = TRUE\n";
*/
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "     -> Basepoint is valid! " << std::endl;
					#endif
					if(  d==dx  ){
						score = basepoint_runner->dynscore - (F*(current->i-basepoint_runner->i));
					}else{
						g = abs(d-dx)-1;
						if( d<dx ) z = current->i - basepoint_runner->i;
						else       z = current->j - basepoint_runner->j;
						score = basepoint_runner->dynscore - (G+(E*g))-(F*z);
					}
					if( score>max_score ){
						max_score = score;
						tmp = basepoint_runner;
					}
				}
			}
			//check if current is in window 
			if( (dx-d)>delta ){ 
				#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "     Breaking here, successing basepoints will not be in the delta-window!";
						std::cerr << std::endl;
				#endif
				break; 
			}
			basepoint_runner = basepoint_runner->inter_d_next;	
		}

		if( max_score==-INFINITY ){
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << "     There is no basepoint in this delta-window!!!" << std::endl;
				std::cerr << "     This message must not appear when grid-points were added." << std::endl;
			#endif
			// no basepoint in window
			max_score = current->kmerscore;
			tmp       = 0;
		}else{
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << "     There is a winner, i:" << tmp->i << " j:" << tmp->j;
				std::cerr << " dynscore:" << max_score;
				
			#endif
			max_score += ( std::min(_k, current->i-tmp->i) / k ) * current->kmerscore;	
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << " dynscore+current-kmer-score:" << max_score << std::endl;
			#endif
		}

		if( current->kmerscore>max_score ){
			max_score = current->kmerscore;
			tmp       = 0;
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << "     It is better to start here, i:" << current->i << " j:" << current->j;
				std::cerr << " dynscore:" << max_score << std::endl;	
			#endif
		}
// -------------------- print -----------------------------	
/*		if (current != 0 && tmp != 0){
			std::cout << "(" << current->i << "," << current->j << ") ---> (" << tmp->i << "," << tmp->j << ")\n";
			std::cout << current << " " << tmp << "\n";
			std::cout << "score = " << max_score << "\n";
			std::cout << "dummy: "<<current->dummy << ", " << tmp->dummy << "\n";
		}
		else{
			if (current == 0) std::cout << "current = 0\n";
			else std::cout << "current = " << "(" << current->i << "," << current->j << ")\n";
			if (tmp == 0) std::cout << "tmp = 0\n";
			else std::cout << "tmp = " << tmp->i << "," << tmp->j << ")\n";
		}
		std::cout << "\n";*/
// ----------------- print ende ---------------------------		
		current->dynscore  = max_score;	
		current->traceback = tmp;

		if( overall_max_score<max_score ){
			overall_max_score = max_score;
			_trace             = current;	
		}
			
		//update head
		#ifdef __CM_DEBUG_KMERALN_DYNAMIC
			std::cerr << "   Updating the head of the trailing list which keeps all matches within the last " << delta;
			std::cerr << " positions (with respect to 'i') ..." << std::endl;
		#endif
		ihead = current;
		//set basepoint_runner to the lowest basepoint in the segment
		basepoint_runner = diag_ptr[segment];
		if( basepoint_runner==0 ){
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << "     Current("<< current <<") is the first basepoint for segment:" << segment << std::endl;
			#endif
			current->inter_d_next = 0;
			current->inter_d_prev = 0;
			current->intra_d_next = 0;
			current->intra_d_prev = 0;
			segmentx = segment;
			while( segmentx >=0 ){
				#ifdef __CM_DEBUG_KMERALN_DYNAMIC
					std::cerr << "          Updating segment " << segmentx << " of addressing array." << std::endl;
				#endif
				//fulfilled for at least the first iteration
				//hence, connection in the else case must be correct
				if( diag_ptr[segmentx]==0 ) diag_ptr[segmentx]=current;
				else{
					//move to the last basepoint!!! in this segment and connect the basepoints
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "          Moving to last basepoint of this segment: " << segmentx;
						std::cerr << std::endl;
					#endif
					tmp = diag_ptr[segmentx ];
					while( tmp->inter_d_next!=0 ) tmp = tmp->inter_d_next;
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "          Connecting " << tmp << " with " << current ;
						std::cerr << " and vice versa" << std::endl;
					#endif
// !!!
					tmp->inter_d_next     = current;
// !!!
					current->inter_d_prev = tmp;
					break;
				}
				--segmentx ;
			}
		}else{
			//find basepoint on the same diagonal or the very next basepoint with a higher diagonal
			//or zero if current is a new basepoint at the end of this segment and no basepoints in upper segments exist
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr <<"     Diagonal+offset of 'current':" << d+offset << std::endl;
			#endif
			while( basepoint_runner!=0 ){
				dx = basepoint_runner->i-basepoint_runner->j;
				#ifdef __CM_DEBUG_KMERALN_DYNAMIC
					std::cerr <<"     Diagonal+offset of basepoint(" << basepoint_runner;
					std::cerr << "):" << dx+offset << std::endl;
				#endif
				if( dx >= d ) break;
				tmp              = basepoint_runner;	
				basepoint_runner = basepoint_runner->inter_d_next;	
			}
			if( basepoint_runner==0 ){
				//current is a new basepoint at the end of the segment
				//diag_ptr is up to date
				#ifdef __CM_DEBUG_KMERALN_DYNAMIC
					std::cerr << "     New basepoint in segment " << segment;
					std::cerr << " without successors in higher segments!" <<std::endl;
				#endif
				tmp->inter_d_next     = current;
				current->inter_d_next = 0;
				current->inter_d_prev = tmp;
				current->intra_d_next = 0;
				current->intra_d_prev = 0;
			}else{
				if( dx==d ){
					//current replaces a basepoint
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "     'current' is on a diagonal with other matches, segment:";
						std::cerr << segment << std::endl;
						std::cerr << "     Updating basepoint:" << basepoint_runner << std::endl;
					#endif
					current->intra_d_prev = 0;	
					current->inter_d_next = basepoint_runner->inter_d_next;
					current->inter_d_prev = basepoint_runner->inter_d_prev;
					if( basepoint_runner->inter_d_next!=0 ){
						basepoint_runner->inter_d_next->inter_d_prev = current;
					}
					if( basepoint_runner->inter_d_prev!=0 ){
						basepoint_runner->inter_d_prev->inter_d_next = current;
					}
					basepoint_runner->intra_d_prev = current;
					current->intra_d_next          = basepoint_runner;
					segmentx                       = segment;
					while( segmentx >=0 && diag_ptr[segmentx ]==basepoint_runner ){
						#ifdef __CM_DEBUG_KMERALN_DYNAMIC
							std::cerr << "          Updating segment " << segmentx;
							std::cerr << " of addressing array." << std::endl;
						#endif
						diag_ptr[segmentx]=current;
						--segmentx;
					}
				}else{
					//current is a new basepoint 
					//insertion before basepoint_runner
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "     Current is inserted in front of the basepoint:";
						std::cerr << basepoint_runner << std::endl;
					#endif
					current->inter_d_prev = basepoint_runner->inter_d_prev;
					if( basepoint_runner->inter_d_prev!=0 ){
						basepoint_runner->inter_d_prev->inter_d_next = current;
					}
					current->inter_d_next          = basepoint_runner;
					basepoint_runner->inter_d_prev = current;
					current->intra_d_next          = 0;
					current->intra_d_prev          = 0;
					if( d<(diag_ptr[segment]->i-diag_ptr[segment]->j) ){
						tmp       = diag_ptr[segment];
						segmentx  = segment;
						while( segmentx >=0 && diag_ptr[segmentx]==tmp ){
							#ifdef __CM_DEBUG_KMERALN_DYNAMIC
								std::cerr << "          Updating segment " << segmentx;
								std::cerr << " of addressing array." << std::endl;
							#endif
							diag_ptr[segmentx]=current;
							--segmentx;
						}
					}
				}
			}
		}
		//update tail
		#ifdef __CM_DEBUG_KMERALN_DYNAMIC
			std::cerr << "   Updating the tail of the trailing list, some matches may fall out ... " << std::endl;
		#endif
		while( itail!=ihead ){
			#ifdef __CM_DEBUG_KMERALN_DYNAMIC
				std::cerr << "      Checking " << itail << " i:" <<  itail->i << " j:" <<  itail->j << std::endl;
			#endif
			if( itail->i < (current->i-delta) ){
				#ifdef __CM_DEBUG_KMERALN_DYNAMIC
					std::cerr << "        --> Remove it:" << std::endl;
				#endif
				if( itail->intra_d_prev==0 ){
					//_itail is a basepoint
					if( itail->inter_d_prev!=0 ){
						itail->inter_d_prev->inter_d_next = itail->inter_d_next;
					}
					if( itail->inter_d_next!=0 ){
						itail->inter_d_next->inter_d_prev = itail->inter_d_prev;
					}
					//check diagonal pointers
					dx       = itail->i - itail->j;
					segmentx = (dx+offset)/W;
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "        itail:" << itail << " is a basepoint in segment ";
						std::cerr << segmentx << std::endl;
						std::cerr << "        Connecting previous(" << itail->inter_d_prev;
						std::cerr << ") and next(" << itail->inter_d_next << ")" << std::endl;
					#endif
					while( segmentx>=0 && diag_ptr[segmentx]==itail ){
						#ifdef __CM_DEBUG_KMERALN_DYNAMIC
							std::cerr << "          Updating segment " << segmentx;
							std::cerr << " of addressing array." << std::endl;
						#endif
						diag_ptr[segmentx]=itail->inter_d_next;
						--segmentx ;
					}
				}else{
					//itail is somewhere on a diagonal with other matches
					#ifdef __CM_DEBUG_KMERALN_DYNAMIC
						std::cerr << "        itail:" << itail << " is on diagonal ";
						std::cerr << itail->i-itail->j+offset;
						std::cerr << " ->setting intra_d_next pointer of previous to zero." << std::endl;
					#endif
					itail->intra_d_prev->intra_d_next=0;
				}
			}else{
				break;	
			}	
			itail=itail->inext;	
		}
		current = current->inext;
	}

	#ifdef __CM_DEBUG_KMERALN_DYNAMIC
		std::cerr << " -----------------------------------------------------------------------------------------" << std::endl;
		std::cerr << " Best cell is i=" << _trace->i << " j=" << _trace->j << " with score= " << _trace->dynscore << std::endl;
		std::cerr << " Number of visited basepoints            : " << visited_basepoints  << std::endl;
		std::cerr << " Number of visited inter-diagonal points : " << visited_interpoints << std::endl;
	#endif
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "Dynamic programming on kmer-match-list done." << std::endl;
	#endif
}


void Kmeraln::_fill(){
	const int offset = Ly-1;

	Kmeraln::_extended_match *current          = _first;
	Kmeraln::_extended_match *itail            = _first;
	Kmeraln::_extended_match *ihead            = 0;
	Kmeraln::_extended_match *basepoint_runner = 0;
	Kmeraln::_extended_match *intra_runner     = 0;
	Kmeraln::_extended_match *tmp              = 0;

	const int mid = (Lx+Ly-1)/2;
	Kmeraln::_extended_match *low_d_basepoint  = 0;
	Kmeraln::_extended_match *high_d_basepoint = 0;
	bool forward = true;
	//inter-diagonal gap length
	int g;
		
	//intra-diagonal (pseudo) gap length
	int z;

	const float k = _k;

	register int d,dx=0;
	float max_score, score;
	
	
	float overall_max_score = -INFINITY;

	while( current!=0 ){

		d      = current->i - current->j;

		//decide where to start
		if( d+offset > mid ) forward=false;
		if(forward)
			basepoint_runner = low_d_basepoint;
		else
			basepoint_runner = high_d_basepoint;
		max_score = -INFINITY;
		while( basepoint_runner!=0 ){
			dx = basepoint_runner->i-basepoint_runner->j;
			if( basepoint_runner->j>current->j ){
				intra_runner = basepoint_runner->intra_d_next;
				while( intra_runner!=0 ){
					if( _in_window(current, intra_runner) ){
						dx = intra_runner->i-intra_runner->j;
						g  = abs(d-dx)-1;
						z  = current->j-intra_runner->j;
						score = intra_runner->dynscore - (G+(E*g)) - (F*z);
						if( score>max_score ){
							max_score = score;
							tmp       = intra_runner;
						}
						break;
					}
					intra_runner = intra_runner->intra_d_next;
				}
			}else{
				if( _in_window(current, basepoint_runner) ){
					if(  d==dx  ){
						score = basepoint_runner->dynscore - (F*(current->i-basepoint_runner->i));
					}else{
						g = abs(d-dx)-1;
						if( d<dx ) z = current->i-basepoint_runner->i;
						else       z = current->j-basepoint_runner->j;
						score = basepoint_runner->dynscore - (G+(E*g))-(F*z);
					}
					if( score>max_score ){
						max_score = score;
						tmp = basepoint_runner;
					}
				}
			}
			//check if current is in window 
			if( forward && (dx-d)>delta ){ break; }
			if( !forward && (d-dx)>delta ){ break; }
			if(forward)
				basepoint_runner = basepoint_runner->inter_d_next;
			else	
				basepoint_runner = basepoint_runner->inter_d_prev;
		}

		if( max_score==-INFINITY ){
			// no basepoint in window
			max_score = current->kmerscore;
			tmp       = 0;
		}else{
			max_score += ( std::min(_k, current->i-tmp->i) / k ) * current->kmerscore;	
		}

		if( current->kmerscore>max_score ){
			max_score = current->kmerscore;
			tmp       = 0;
		}

		current->dynscore  = max_score;	
		current->traceback = tmp;
		if( overall_max_score<max_score ){
			overall_max_score = max_score;
			_trace             = current;	
		}

		ihead = current;
		if(forward)
			basepoint_runner = low_d_basepoint;
		else
			basepoint_runner = high_d_basepoint;
		if( basepoint_runner==0 ){
			current->inter_d_next = 0;
			current->inter_d_prev = 0;
			current->intra_d_next = 0;
			current->intra_d_prev = 0;
			low_d_basepoint       = current;
			high_d_basepoint      = current;
		}else{
			if(forward){
				while( basepoint_runner!=0 ){
					dx = basepoint_runner->i-basepoint_runner->j;
					if( dx >= d ) break;
					tmp              = basepoint_runner;	
					basepoint_runner = basepoint_runner->inter_d_next;	
				}
				
				if( basepoint_runner==0 ){
					tmp->inter_d_next     = current;
					current->inter_d_next = 0;
					current->inter_d_prev = tmp;
					current->intra_d_next = 0;
					current->intra_d_prev = 0;
					high_d_basepoint      = current;
				}else{
					if( dx==d ){
						current->intra_d_prev = 0;	
						current->inter_d_next = basepoint_runner->inter_d_next;
						current->inter_d_prev = basepoint_runner->inter_d_prev;
						if( basepoint_runner->inter_d_next!=0 ){
							basepoint_runner->inter_d_next->inter_d_prev = current;
						}
						if( basepoint_runner->inter_d_prev!=0 ){
							basepoint_runner->inter_d_prev->inter_d_next = current;
						}
						basepoint_runner->intra_d_prev = current;
						current->intra_d_next          = basepoint_runner;
						if( low_d_basepoint==basepoint_runner ) 
							low_d_basepoint = current;
						if( high_d_basepoint==basepoint_runner ) 
							high_d_basepoint = current;
					}else{
						current->inter_d_prev = basepoint_runner->inter_d_prev;
						if( basepoint_runner->inter_d_prev!=0 ){
							basepoint_runner->inter_d_prev->inter_d_next = current;
						}else{
							low_d_basepoint = current;
						}
						current->inter_d_next          = basepoint_runner;
						basepoint_runner->inter_d_prev = current;
						current->intra_d_next          = 0;
						current->intra_d_prev          = 0;
					}
				}
			}else{
				while( basepoint_runner!=0 ){
					dx = basepoint_runner->i-basepoint_runner->j;
					if( dx <= d ) break;
					tmp              = basepoint_runner;	
					basepoint_runner = basepoint_runner->inter_d_prev;	
				}
				if( basepoint_runner==0 ){
					tmp->inter_d_prev     = current;
					current->inter_d_next = tmp;
					current->inter_d_prev = 0;
					current->intra_d_next = 0;
					current->intra_d_prev = 0;
					low_d_basepoint      = current;
				}else{
					if( dx==d ){
						current->intra_d_prev = 0;	
						current->inter_d_next = basepoint_runner->inter_d_next;
						current->inter_d_prev = basepoint_runner->inter_d_prev;
						if( basepoint_runner->inter_d_next!=0 ){
							basepoint_runner->inter_d_next->inter_d_prev = current;
						}
						if( basepoint_runner->inter_d_prev!=0 ){
							basepoint_runner->inter_d_prev->inter_d_next = current;
						}
						basepoint_runner->intra_d_prev = current;
						current->intra_d_next          = basepoint_runner;
						if( low_d_basepoint==basepoint_runner ) 
							low_d_basepoint = current;
						if( high_d_basepoint==basepoint_runner ) 
							high_d_basepoint = current;
					}else{
						current->inter_d_prev = basepoint_runner;
						current->inter_d_next = basepoint_runner->inter_d_next;
						if( basepoint_runner->inter_d_next!=0 ){
							basepoint_runner->inter_d_next->inter_d_prev = current;
						}else{
							high_d_basepoint = current;
						}
						basepoint_runner->inter_d_next = current;
						current->intra_d_next          = 0;
						current->intra_d_prev          = 0;
					}
				}
			}
		}	
		//update tail
		while( itail!=ihead ){
			if( itail->i < (current->i-delta) ){
				if( itail->intra_d_prev==0 ){
					//_itail is a basepoint
					if( itail->inter_d_prev!=0 ){
						itail->inter_d_prev->inter_d_next = itail->inter_d_next;
					}
					if( itail->inter_d_next!=0 ){
						itail->inter_d_next->inter_d_prev = itail->inter_d_prev;
					}
					if( itail==low_d_basepoint && itail==high_d_basepoint){
						low_d_basepoint  = itail->inter_d_next;
						high_d_basepoint = low_d_basepoint;
					}else if( itail==high_d_basepoint ){
						high_d_basepoint = high_d_basepoint->inter_d_prev;
					}else if( itail==low_d_basepoint ){
						low_d_basepoint = low_d_basepoint->inter_d_next;
					}
				}else{
					itail->intra_d_prev->intra_d_next=0;
				}
			}else{
				break;	
			}	
			itail=itail->inext;	
		}
		current = current->inext;
	}
}

float Kmeraln::get_score(){
	return _score;
}

float Kmeraln::get_rel_score(){
	return _score/std::min(Lx,Ly);
}

float Kmeraln::get_evalue(){
	float query_len_norm = Lx - (log(K * Lx * Ly) / H);
	float db_len_norm = db_len - (db_seqs * log (K * Lx) / H) - (sum_log_db_len / H);
	
/*	std::cerr << std::endl << "Query length : " << Lx << std::endl;
	std::cerr << "Query length norm. : " << query_len_norm << std::endl;
	std::cerr << "DB length : " << db_len << std::endl;
	std::cerr << "DB length norm. : " << db_len_norm << std::endl;
	std::cerr << "Score : " << _score << std::endl;
	std::cerr << "E-value : " << (0.5 * K * query_len_norm * db_len_norm * exp(-Lamda * _score)) << std::endl;
*/	
	return 0.5 * K * query_len_norm * db_len_norm * exp(-Lamda * _score);
}

float Kmeraln::get_aln_cov_short(){
	return aln_cov_short;
}

float Kmeraln::get_aln_cov_long(){
	return aln_cov_long;
}

size_t Kmeraln::get_number_of_matches(){
	return ExtMatchHeap->get_number_of_current_objects_in_use();
}


inline bool Kmeraln::_in_window(Kmeraln::_extended_match *current, Kmeraln::_extended_match *candidate){
	//allow gridpoints 

/*	if( current->i==candidate->i && current->j==candidate->j ) return true;
	if( current->dummy && candidate->dummy ) return true;
	if( current->j <= candidate->j ) return false;
	if( current->j-delta > candidate->j ) return false;
	if( current->i <= candidate->i ) return false;
	return true;*/

/*    if (current->i == candidate->i && current->j == candidate->j) 
    	return true;
	else if (current->dummy && candidate->dummy)
		return (current->i >= candidate->i) && (current->j >= candidate->j) && (current->i - delta <= candidate->i) && (current->j - delta <= candidate->j);
	else
		return (current->i > candidate->i) && (current->j > candidate->j) && (current->i - delta <= candidate->i) && (current->j - delta <= candidate->j);
*/
	if( current->i==candidate->i && current->j==candidate->j ) return true;
	if( current->dummy && candidate->dummy ){ 
		if( current->j < candidate->j ) return false;
		if( current->j-delta > candidate->j ) return false;
		if( current->i < candidate->i ) return false;
	}
	else{
		if( current->j <= candidate->j ) return false;
		if( current->j-delta > candidate->j ) return false;
		if( current->i <= candidate->i ) return false;
	}
	return true;
}



Pair_aln* Kmeraln::traceback(){
	#ifdef __CM_DEBUG_KMERALN
		std::cerr << "Traceback ..." << std::endl;
	#endif
	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		int traceback_cell_count = 0;
		int dummy_cell_count     = 0;
	#endif
	if( _trace==0 ) return new Pair_aln(1);

	Pair_aln *ret = new Pair_aln(Lx+Ly);
	ret->len_x     = Lx;
	ret->len_y     = Ly;
	int alnpos     = Lx+Ly-1;

	Kmeraln::_extended_match *tmp  = _trace;
	Kmeraln::_extended_match *prev = _trace;
	_score                         = 0.0f;
	int i,j,dc,dp,g,z,istop;
	
	//in the first iteration only the trailing segment of the start traceback kmer is processed
	i = _trace->i + _k -1;
	j = _trace->j + _k -1;
	
	while( tmp != 0 ){
	
		if( tmp->dummy ) {
			tmp = tmp->traceback;
			continue;
		}

		dc     = tmp->i  - tmp->j;
		dp     = prev->i - prev->j;
		

		if (prev->i < tmp->i) {
			std::cout << "\n!!! prev->i = " << prev->i << ", tmp->i = " << tmp->i << "\n";
			std::cout << "    prev->j = " << prev->j << ", tmp->j = " << tmp->j << "\n";
		}
		if (prev->j < tmp->j) {
			std::cout << "\n    prev->i = " << prev->i << ", tmp->i = " << tmp->i << "\n";
			std::cout << "!!! prev->j = " << prev->j << ", tmp->j = " << tmp->j << "\n\n";
		}
		
		if( dp < dc ){
			//gap in x -> insert in y
			z = prev->j - tmp->j;
			g = abs(dp-dc);
			if ((z-g) < 0){
/*			++traceback_cell_count;
			std::cerr << "                       Match  i:" << tmp->i;
			std::cerr << " j:" << tmp->j << " dp:" << tmp->i - tmp->j << "\n";
			std::cerr << "                       Prev   i:" << prev->i;
			std::cerr << " j:" << prev->j << " dc:" << prev->i - prev->j << "\n";
			std::cerr << " score:" << tmp->kmerscore << " dynscore:" << tmp->dynscore ;
			std::cerr << " current aln score:" << _score << std::endl;*/
			}
			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "                       Gap in X  g:" << g << " z:" <<  z << std::endl;
			#endif
			const float tmp_score = _compute_gapx( i, j, z, g, ret, alnpos);
			_score += tmp_score;
			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "                       gap score:" << tmp_score << std::endl;
				std::cerr << "                           score:" << _score << std::endl;
			#endif
			i = tmp->i + _kmer_offset;
			j = tmp->j + _kmer_offset;
			
			
			
			
		}else if( dp>dc ){		
			//gap in y -> insert in x
			z = prev->i - tmp->i;
			g = abs(dp-dc);
			if ((z-g) < 0){
/*			++traceback_cell_count;
			std::cerr << "                       Match  i:" << tmp->i;
			std::cerr << " j:" << tmp->j << " dp:" << tmp->i - tmp->j << "\n";
			std::cerr << "                       Prev   i:" << prev->i;
			std::cerr << " j:" << prev->j << " dc:" << prev->i - prev->j << "\n";
			std::cerr << " score:" << tmp->kmerscore << " dynscore:" << tmp->dynscore ;
			std::cerr << " current aln score:" << _score << std::endl;*/
			}
			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "                       Gap in Y  g:" << g << " z:" <<  z << std::endl;
			#endif
			const float tmp_score = _compute_gapy(i, j, z, g, ret, alnpos);
			_score += tmp_score;
/*			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "                       gap score:" << tmp_score << std::endl;
				std::cerr << "                           score:" << _score << std::endl;
			#endif*/
			i = tmp->i +  _kmer_offset;
			j = tmp->j +  _kmer_offset;
			
			
			
			
		}else if( dp ==dc ){
			istop = tmp->i+_kmer_offset;
			while( i>istop ){
//				#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
//					std::cerr << "   Aligning " <<  int2aa[X[i]] << "(" << i <<")   " <<  int2aa[Y[j]] << "(" << j << ")"<< std::endl;
//				#endif
//				_score += matrix[ X[i] ][ Y[j] ];
				_score += X_seq->get_score_with(i, Y[j]);
				ret->x[alnpos] = i;
				ret->y[alnpos] = j;
				--i;--j;--alnpos;
			}
		}
		prev = tmp;
		tmp = tmp->traceback;
	}
	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << "  -----------------------------------------------------------------" << std::endl;
		std::cerr << "  Raw alignment score  without extension: " << _score << std::endl;
	#endif

	//extend ends
	//move from the beginning of the alignment towards the left side
	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << "  Extension on diagonal to the left side:" << std::endl;
		std::cerr << "     Start at i:" << i << " j:" << j << "     score:" << _score << std::endl;
	#endif
	int save_i              = i;
	int save_j              = j;
	float cur_run_out_score = _score;
	float max_run_out_score = _score;
	int   max_run_out_pos   = i+1;
	while( i>=0 && j>=0 ){
//		cur_run_out_score += matrix[ X[i] ][ Y[j] ];
		cur_run_out_score += X_seq->get_score_with(i, Y[j]);
		if( cur_run_out_score>max_run_out_score  ){ 
			max_run_out_score = cur_run_out_score ;
			max_run_out_pos   = i;
		}else if( (max_run_out_score-cur_run_out_score)>=_dropoff ) break;
		--i; --j;
	}

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << "      Maximum score reached at:" << max_run_out_pos << "  score:" << max_run_out_score << std::endl;
	#endif

	i                    = save_i;
	j                    = save_j;
	while( i>=max_run_out_pos ){
		ret->x[alnpos] = i;
		ret->y[alnpos] = j;
		--i; --j;--alnpos;
	}

	ret->x_aln_start = i+2;
	ret->y_aln_start = j+2;

	//move alignment to the left in the arrays
	int alnlen = Lx+Ly-alnpos-1;
	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << std::endl;
		std::cerr << "  Moving alignment data in array to the left, length of alignment is " << alnlen << std::endl;
		std::cerr << std::endl;
	#endif
	if( alnlen!=0 ){
		memmove(ret->x, &ret->x[alnpos+1], sizeof(int)*alnlen);
		memmove(ret->y, &ret->y[alnpos+1], sizeof(int)*alnlen);
	}

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << "  Extension on diagonal to the right:" << std::endl;
		std::cerr << "     Start at i:" << _trace->i + _k << " j:" << _trace->j + _k << "     score:" << _score << std::endl;
	#endif

	//extend to the right side
	i                   = _trace->i + _k;
	j                   = _trace->j + _k;
	max_run_out_score   = _score;
	cur_run_out_score   = _score;
	max_run_out_pos     = i-1;
	while( i<Lx && j<Ly ){
//		cur_run_out_score += matrix[ X[i] ][ Y[j] ];
		cur_run_out_score += X_seq->get_score_with(i, Y[j]);
		if( cur_run_out_score>max_run_out_score ){
			max_run_out_score = cur_run_out_score;
			max_run_out_pos   = i;
		}else if( (max_run_out_score-cur_run_out_score)>=_dropoff ) break;
		++i; ++j;
	}

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << "      Maximum score reached at:" << max_run_out_pos << "  score:" << max_run_out_score << std::endl;
	#endif

	i = _trace->i + _k;
	j = _trace->j + _k;
	while( i<=max_run_out_pos ){
		ret->x[alnlen] = i;
		ret->y[alnlen] = j;
		++i; ++j;++alnlen;
	}

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << std::endl;
		std::cerr << "   Length of alignment:" << alnlen << std::endl;
		std::cerr << std::endl;
	#endif

	ret->x_aln_end = i;
	ret->y_aln_end = j;

	//compute best aln-segment
	//left_pos[i] holds the position before the next positive segment
	int right_pos_max      = 0;
	int left_pos_max       = 0;
	float right_score_max  = 0.0f;	
	float tmp_score        = 0.0f;
	int left_start         = -1;
	bool gap               = false;

	for( int li=0; li<alnlen; ++li){
		if( ret->x[li]==-1 ){
			if( gap ){
				tmp_score -= E;
			}else{
				tmp_score -= G;
				gap = true;	
			}
		}else if( ret->y[li]==-1 ){
			if( gap ){
				tmp_score -= E;
			}else{
				tmp_score -= G;
				gap = true;	
			}
		}else{
//			tmp_score += matrix[ X[ret->x[li]] ][ Y[ret->y[li]] ];
			tmp_score += X_seq->get_score_with(ret->x[li], Y[ret->y[li]]);
			gap = false;
		}
	
		if( tmp_score < 0.0f ){
			left_start = li;
			tmp_score = 0.0f;
		}
		if( right_score_max<tmp_score ){
			right_score_max = tmp_score;
			right_pos_max   = li;
			left_pos_max    = left_start+1;
		}
/*		#ifdef __CM_DEBUG_KMERALN_TRACEBACK
			char xc = '-';
			char yc = '-';
			if( ret->x[li]!=-1 ) xc = int2aa[ X[ret->x[li]] ];
			if( ret->y[li]!=-1 ) yc = int2aa[ Y[ret->y[li]] ];
			std::cerr << li << " " << tmp_score  << "   " << xc << " " << yc << std::endl;
		#endif*/
	}
	
	if( right_pos_max!=0 ) alnlen = right_pos_max-left_pos_max+1;

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK
		std::cerr << "   Starting position at the left side  : " << left_pos_max << std::endl;
		std::cerr << "        Belongs to position (in x)     : " << ret->x[left_pos_max] << std::endl;
		std::cerr << "        Belongs to position (in y)     : " << ret->y[left_pos_max] << std::endl;
		std::cerr << std::endl;
		std::cerr << "   Starting position at the right side : " << right_pos_max << std::endl;
		std::cerr << "         Belongs to position (in x)    : " << ret->x[right_pos_max] << std::endl;
		std::cerr << "         Belongs to position (in y)    : " << ret->y[right_pos_max] << std::endl;
		std::cerr << std::endl;
		std::cerr << "   Score for alignment                 : " << right_score_max << std::endl;
		std::cerr << "   Length of alignment                 : " << alnlen << std::endl;
	#endif

	if( alnlen!=0 ){
		memmove(ret->x, &ret->x[left_pos_max], sizeof(int)*alnlen);
		memmove(ret->y, &ret->y[left_pos_max], sizeof(int)*alnlen);
	}

	_score                 = right_score_max;
	ret->len_aln           = alnlen;
	ret->score             = _score;
	for(int t=0; t<alnlen; ++t){
		if( ret->x[t]==-1 ){
			++ret->gaps;
			++ret->gaps_x;
			++ret->len_equ;
//			++ret->aln_len_y;
		}else if( ret->y[t]==-1 ){
			++ret->gaps;
			++ret->gaps_y;
			++ret->len_equ;
//			++ret->aln_len_x;
		}else{ 
			char xc = X[ret->x[t]];
			char yc = Y[ret->y[t]];
			if( xc==yc && int2aa[xc]!='X' ) {
				++ret->idents;
				++ret->pos_matches;
			}
//			else if( matrix[ xc ][ yc ]>0.0f ) ++ret->pos_matches;
			else if (X_seq->get_score_with(ret->x[t], yc) > 0.0f) ++ret->pos_matches;
			else ++ret->neg_zero_matches;
			++ret->len_equ;
			++ret->aln_len_x;
			++ret->aln_len_y;
		}
	}
	
	if( Lx<Ly ){
		aln_cov_short = ret->aln_len_x/((float)Lx + (float)ret->gaps_x);
		aln_cov_long  = ret->aln_len_y/((float)Ly + (float)ret->gaps_y);
	}else{
		aln_cov_short = ret->aln_len_y/((float)Ly + (float)ret->gaps_y);
		aln_cov_long  = ret->aln_len_x/((float)Lx + (float)ret->gaps_x);
	}

	#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
		std::cerr << " ----------------------------------------" << std::endl;
		std::cerr << " Number of traceback cells : " << traceback_cell_count << std::endl;
		std::cerr << " Number of dummy cells     : " << dummy_cell_count     << std::endl;
	#endif
	#ifdef __CM_DEBUG_KMERALN 
		std::cerr << "Traceback done." << std::endl;
	#endif

	return ret;
}

inline float Kmeraln::_compute_gapx(	const int i, 
						const int j, 
						const int z, 
						const int g,
						Pair_aln *alnp, 
						int &alnpos){
	float ret   = -INFINITY;
	int places  = z-g;
	if (places < 0) std::cout << "(X) Places: " << places << "\n";
	int i_f     = i - places + 1;
	int j_f     = j - z + 1;
	int i_b     = i;
	int j_b     = j;
	float tmp_f = 0.0f;
	float tmp_b = 0.0f;

	forward[0] = 0.0f;
	backward[places] = 0.0f;

	for( int t=0; t<places; ++t){
//		tmp_f += matrix[ X[i_f] ][ Y[j_f] ];
//		tmp_b += matrix[ X[i_b] ][ Y[j_b] ];
		tmp_f += X_seq->get_score_with(i_f, Y[j_f]);
		tmp_b += X_seq->get_score_with(i_b, Y[j_b]);
		forward[t+1]         = tmp_f;
		backward[places-t-1] = tmp_b; 
		++i_f; ++j_f; --i_b; --j_b;
	}

	int tmp_t = 0;
	for( int t=0; t<=places; ++t){
		if( forward[t] + backward[t] > ret ){
			ret   = forward[t] + backward[t];
			tmp_t = t;
		}
	}
	ret -= (G+(E*(g-1)));

	i_b     = i;
	j_b     = j;
	for( int t=places; t>=0; --t){
		if( t == tmp_t ){
				for( int gi=0; gi<g; ++gi){
/*					#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
//						std::cerr << "   Aligning " <<  "-" << "   " <<  int2aa[Y[j_b]] << std::endl;
					#endif*/
					alnp->x[alnpos] = -1;
					alnp->y[alnpos] = j_b;
					--alnpos;
					--j_b;
				}
		}else{
/*			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "   Aligning " <<  int2aa[X[i_b]] << "   " <<  int2aa[Y[j_b]] << std::endl;
			#endif*/
			alnp->x[alnpos] = i_b;
			alnp->y[alnpos] = j_b;
			--alnpos;
			--i_b;
			--j_b;
		}
	}
	return ret;
}

inline float Kmeraln::_compute_gapy(	const int i, 
						const int j, 
						const int z, 
						const int g, 
						Pair_aln *alnp, 
						int &alnpos){
	float ret   = -INFINITY;
	int places  = z-g;
	if (places < 0) std::cout << "(Y) Places: " << places << "\n";
	int i_f     = i - z + 1;
	int j_f     = j - places + 1;
	int i_b     = i;
	int j_b     = j;
	float tmp_f = 0.0f;
	float tmp_b = 0.0f;

	forward[0] = 0.0f;
	backward[places] = 0.0f;

	for( int t=0; t<places; ++t){
//		tmp_f += matrix[ X[i_f] ][ Y[j_f] ];
//		tmp_b += matrix[ X[i_b] ][ Y[j_b] ];
		tmp_f += X_seq->get_score_with(i_f, Y[j_f]);
		tmp_b += X_seq->get_score_with(i_b, Y[j_b]);
		forward[t+1]         = tmp_f;
		backward[places-t-1] = tmp_b; 
		++i_f; ++j_f; --i_b; --j_b;
	}
	
	int tmp_t = 0;
	for( int t=0; t<=places; ++t){
		if( forward[t] + backward[t] > ret ){
			ret   = forward[t] + backward[t];
			tmp_t = t;
		}
	}
	ret -= (G+(E*(g-1)));

	i_b     = i;
	j_b     = j;
	for( int t=places; t>=0; --t){
		if( t == tmp_t ){
				for( int gi=0; gi<g; ++gi){
/*					#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
						std::cerr << "   Aligning " <<  int2aa[X[i_b]] << "   " <<  "-" << std::endl;
					#endif*/
					alnp->x[alnpos] = i_b;
					alnp->y[alnpos] = -1;
					--alnpos;
					--i_b;
				}
		}else{
/*			#ifdef __CM_DEBUG_KMERALN_TRACEBACK 
				std::cerr << "   Aligning " <<  int2aa[X[i_b]] << "   " <<  int2aa[Y[j_b]] << std::endl;
			#endif*/
			alnp->x[alnpos] = i_b;
			alnp->y[alnpos] = j_b;
			--alnpos;
			--i_b;
			--j_b;
		}
	}
	return ret;
}

std::ostream &Kmeraln::print_trace(std::ostream &out){
	Kmeraln::_extended_match *tmp = _trace;
	while( tmp!=0 ){
		std::cout <<  tmp->i << " " << tmp->j << std::endl;
		tmp = tmp->traceback;
	}
	return out;
}

std::ostream &Kmeraln::print_vertex_list(std::ostream &out){
	Kmeraln::_extended_match *tmp=_first;
	uint pos = 0;
	while( tmp!=0 ){
		if( tmp->dummy)
			std::cout << tmp->i << " " << tmp->j << " " << tmp->kmerscore << std::endl;
		else
			std::cerr<< tmp->i << " " << tmp->j << " " << tmp->kmerscore << std::endl;
		tmp = tmp->inext;
		++pos;
	}
	return out;
}

/*float Kmeraln::_compute_overlap_score(char *alnx, char *alny, size_t len, const float **matrix){
	
	float ret=0;
	
	int idents    = 0; 
	int gaps      = 0;
	int positives = 0;
	int neg_zeros = 0;
	int len_aln   = 0;
	int len_equ   = 0;

	size_t start=0;
	for(size_t i=0; i<len; ++i){
		if( (alnx[i]=='-') || (alny[i]=='-') ) {
			++start;
			++len_aln;
		}else break;
	}
	size_t stop=len;
	for(size_t i=len-1; i>=0; --i){
		if( (alnx[i]=='-') || (alny[i]=='-') ){ 
			--stop;
			++len_aln;
		}else break;
	}
	//std::cout << start << " " <<stop << std::endl;
	bool open = false;
	for(size_t i=start; i<stop; ++i){
		//std::cout << "RET:" << ret << std::endl;
		if(alnx[i]=='-' || alny[i]=='-'){
			++gaps;
			++len_equ;
			if(!open) ret -= G;
			else ret -= E;
			open = true;
		}else{ 
			//std::cout << matrix[ alnx[i] ][ alny[i] ] << std::endl;
			ret += matrix[ alnx[i] ][ alny[i] ];
			if( alnx[i]== alny[i] ) ++idents;
			else if( matrix[ alnx[i] ][ alny[i] ]>0.0f ) ++positives;
			else ++neg_zeros;
			++len_equ;
			open=false;
		}
	}
	len_aln += len_equ;
	std::cerr << "Length of alignment                : " << len_aln   << std::endl;
	std::cerr << "Equivalenced residues (incl. gaps) : " << len_equ   << std::endl;
	std::cerr << "Score (for equiv. residues)        : " << ret       << std::endl;
	std::cerr << "Identities                         : " << idents    << std::endl;
	std::cerr << "Positive matches (incl. idents)    : " << positives << std::endl;
	std::cerr << "Less or equal to zero matches      : " << neg_zeros << std::endl;
	std::cerr << "Gaps                               : " << gaps      << std::endl;

	return ret;
}*/

void Kmeraln::set_coverage_criterion(float coverage){
	COV = coverage-0.001;
}

Pair_aln* Kmeraln::get_aln(){
	return aln;
}
