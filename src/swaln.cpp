/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "swaln.h"


SWaln::SWaln(
	METHOD m, 
	Sequence *sx, 
	Sequence *sy, 
	Matrix *mat,
	const float gapopen, 
	const float gapextend, 
	const int diagonal,
	const unsigned int bandwidth, 
	const bool v): 	
		method(m), 
		L_x(sx->length()+1), 
		L_y(sy->length()+1), 
		x(sx->get_sequence()), 
		y(sy->get_sequence()), 
		AMINOACID_DIM( mat->get_aa_dim() ),
		int2aa( mat->get_int2aa() ),
		aa2int( mat->get_aa2int() ),
		matrix( mat->get_matrix() ),
		G(gapopen), 
		E(gapextend), 
		D(diagonal),
		W(bandwidth),
		verbose(v){
	switch(method){
		case SWaln::OVERLAP_FLIP: _init_overlap_matrices_flip();        break;
		case SWaln::OVERLAP_FULL: _init_overlap_matrices_full();        break;
		case SWaln::OVERLAP_BAND: _init_overlap_matrices_band_flip();   break;
		case SWaln::LOCAL_FLIP:   _init_overlap_matrices_flip();        break;
		case SWaln::LOCAL_BAND:   _init_overlap_matrices_band_flip();   break;
	}
	max_i      = 0;
	max_j      = 0;
	max_j_band = 0;
	score      = 0.0f;
}


SWaln::~SWaln(){
	switch(method){
		case SWaln::OVERLAP_FLIP:{
			_delete_overlap_matrices_flip();
			break;
		}
		case SWaln::OVERLAP_FULL:{
			_delete_overlap_matrices_full();
			break;
		}
		case SWaln::OVERLAP_BAND:{
			_delete_overlap_matrices_band_flip();
			break;
		}
		case SWaln::LOCAL_FLIP:{
			_delete_overlap_matrices_flip();
			break;
		}
		case SWaln::LOCAL_BAND:{
			_delete_overlap_matrices_band_flip();
			break;
		}
	}
}

void SWaln::align(){
	switch(method){
		case SWaln::OVERLAP_FLIP:{
			_align_overlap_flip();
			break;
		}
		case SWaln::OVERLAP_FULL:{
			_align_overlap_full();
			break;
		}
		case SWaln::OVERLAP_BAND:{
			_align_overlap_band_flip( D, W);
			break;
		}
		case SWaln::LOCAL_FLIP:{
			_align_local_flip();
			break;
		}
		case SWaln::LOCAL_BAND:{
			_align_local_band_flip( D, W);
			break;
		}
	}
}


void SWaln::_align_overlap_full() throw (std::exception){
	if(verbose) std::cout << "Fill matrices (OVERLAP_FULL) ..." << std::endl;
	float max_value_i = -INFINITY;
	register int i,j;
	register float m, tx, ty, tmp;
	for( i=1; i<L_x; ++i){
		for( j=1; j<L_y; ++j){
			tmp = matrix[x[i-1]][y[j-1]];
			m   = M[i-1][j-1]   + tmp;
			tx  = I_x[i-1][j-1] + tmp;
			ty  = I_y[i-1][j-1] + tmp;
			if( (m>=tx) && (m>=ty) ){
				M[i][j] = m;
				T[i][j] = 1;
			}else if( (tx>=m) && (tx>=ty) ){
				M[i][j] = tx;
				T[i][j] = 0;
			}else{
				M[i][j] = ty;
				T[i][j] = 2;
			}

			m  = M[i-1][j]   - G;
			tx = I_x[i-1][j] - E;		
			if( m>=tx ){ 
				I_x[i][j] = m;	
			}else{
				I_x[i][j] = tx;
				T[i][j]  += 4;
			}

			m  = M[i][j-1]   - G;
			ty = I_y[i][j-1] - E;		
			if( m>=ty ){ 
				I_y[i][j] = m;		
			}else{
				I_y[i][j] = ty;
				T[i][j]  += 8;
			}
		}
		if(M[i][L_y-1]>max_value_i){
			max_value_i = M[i][L_y-1];
			max_i       = i;
		}
	}
	
	float max_value_j = -INFINITY;
	for( j=0; j<L_y; ++j ){
		if( M[L_x-1][j] > max_value_j ){
			max_value_j = M[L_x-1][j];
			max_j       = j;
		}
	}

	if(max_value_i>max_value_j){
		max_j = L_y-1;
		score = max_value_i;
	}else{
		max_i = L_x-1;
		score = max_value_j;
	}
	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_align_overlap_flip() throw (std::exception){
	if(verbose) std::cout << "Fill matrices (OVERLAP_FLIP) ..." << std::endl;
	float max_value_i = -INFINITY;
	register int i,j;
	register float m, tx, ty, tmp;
	register bool flip=0;
	register const float *MTRX;
	register float *MM, *MM_last;
	register float *II_x, *II_x_last;
	register float *II_y, *II_y_last;
	for( i=1; i<L_x; ++i){
		MM_last    = M[flip];
		II_x_last  = I_x[flip];
		II_y_last  = I_y[flip];
		flip       = !flip;
		MM         = M[flip];
		II_x       = I_x[flip];
		II_y       = I_y[flip];	
		MTRX = matrix[x[i-1]];
		for( j=1; j<L_y; ++j){
			tmp = MTRX[y[j-1]];
			m   = MM_last[j-1]    + tmp;
			tx  = II_x_last[j-1]  + tmp;
			ty  = II_y_last[j-1]  + tmp;		
			if( (m>=tx) && (m>=ty) ){	
				MM[j]   = m;
				T[i][j] = 1;
			}else if( (tx>=m) && (tx>=ty) ){
				MM[j]   = tx;
				T[i][j] = 0;
			}else{
				MM[j]   = ty;
				T[i][j] = 2;
			}

			m  = MM_last[j]   - G;
			tx = II_x_last[j] - E;			
			if( m>=tx ){
				II_x[j] = m;
				
			}else{
				II_x[j] = tx;
				T[i][j] += 4;
			}

			m  = MM[j-1]   - G;
			ty = II_y[j-1] - E;		
			if( m>=ty ){
				II_y[j] = m;
				
			}else{
				II_y[j] = ty;
				T[i][j] += 8;
			}
		}
		
		if(MM[L_y-1]>max_value_i){
			max_value_i = MM[L_y-1];
			max_i       = i;
		}
	}

	float max_value_j = -INFINITY;
	for( j=0; j<L_y; ++j ){
		if( M[flip][j] > max_value_j ){
			max_value_j = M[flip][j];
			max_j       = j;
		}
	}

	if(max_value_i>max_value_j){
		max_j = L_y-1;
		score = max_value_i;
	}else{
		max_i = L_x-1;
		score = max_value_j;
	}
	if(verbose) std::cout << "done." << std::endl;
}

void SWaln::_align_local_flip() throw (std::exception){

	/*__m128 f1,f2;
		
	union _4float{
		__m128 m128;
		float f[4];
	} rf;	*/


	if(verbose) std::cout << "Fill matrices (LOCAL_FLIP) ..." << std::endl;
	score = -INFINITY;
	register int i,j,jj;
	register float m, tx, ty;
	register bool flip=0;

	register const float *MTRX;
	register float *MM, *MM_last;
	register float *II_x, *II_x_last;
	register float *II_y, *II_y_last;
	register float tmp;
	for( i=1; i<L_x; ++i){	
		MM_last    = M[flip];
		II_x_last  = I_x[flip];
		II_y_last  = I_y[flip];
		flip       = !flip;
		MM         = M[flip];
		II_x       = I_x[flip];
		II_y       = I_y[flip];	
		MTRX = matrix[x[i-1]];
		for( j=1; j<L_y; ++j){
			jj  = j-1;
			tmp = MTRX[y[jj]];
			m   = MM_last[jj]    + tmp;
			tx  = II_x_last[jj]  + tmp;
			ty  = II_y_last[jj]  + tmp;
			/*
			f1 = _mm_set_ps( MM_last[jj], II_x_last[jj], II_y_last[jj], 0.0f);
			f2 = _mm_set_ps( tmp, tmp, tmp, 0.0f);
			rf.m128 = _mm_add_ps(f1, f2);
			m  = rf.f[0];
			tx = rf.f[1];
			ty = rf.f[2];*/
			
			if( (m>0) && (m>=tx) && (m>=ty) ){	
				MM[j]   = m;
				T[i][j] = 1;
				if( m>score ){
					score = m;
					max_i = i;
					max_j = j;
				}
			}else if(  (tx>0) && (tx>=m) && (tx>=ty) ){
				MM[j]   = tx;
				T[i][j] = 0;
				if( tx>score ){
					score = tx;
					max_i = i;
					max_j = j;
				}
			}else if(  (ty>0) && (ty>=m) && (ty>=tx) ){
				MM[j]   = ty;
				T[i][j] = 2;
				if( ty>score ){
					score = ty;
					max_i = i;
					max_j = j;
				}
			}else{
				if( 0>score ){
					score = 0;
					max_i = i;
					max_j = j;
				}
				MM[j]   = 0;
				T[i][j] = 3;
			}


			m  = MM_last[j]   - G;
			tx = II_x_last[j] - E;			
			if( m>=tx ){
				II_x[j] = m;
				
			}else{
				II_x[j] = tx;
				T[i][j] += 4;
			}

			m  = MM[jj]   - G;
			ty = II_y[jj] - E;		
			if( m>=ty ){
				II_y[j] = m;
				
			}else{
				II_y[j] = ty;
				T[i][j] += 8;
			}

		}
	}
	//printf("maxi:%i maxj:%i score:%f\n", max_i, max_j, score);
	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_align_overlap_band_flip(const int diagonal, const int bandwidth) throw (std::exception){
	if(verbose) std::cout << "Fill matrices (OVERLAP_BAND) ..." << std::endl;
	if( (diagonal<-((int)L_y-2)) || (diagonal>((int)L_x-2)) ){
		throw MyException("Diagonal(%i) is out of range! Possible values are [%i..%i].\n",diagonal ,-(L_y-2), (L_x-2));
	}	
	score = -INFINITY;
	//main diagonal is 0   d=i-j!
	const int upper_diagonal = diagonal + bandwidth;
	const int lower_diagonal = diagonal - bandwidth;

	register const float *MTRX;
	register float *MM, *MM_last;
	register float *II_x, *II_x_last;
	register float *II_y, *II_y_last;
	// matrix positions, -- i-1 position int the sequence
	// j position in the band
	// j_start start of band idx
	// j_stop last position of band
	register int i, j, jj, j_start, j_stop, j_offset, diag;
	//tmp holds different temporary values e.g. BLOSUM score of amino acids
	//m: BLOSUM score + last match score
	//tx: BLOSUM score + last score of a gap in x 
	//ty: BLOSUM score + last score of a gap in y
	register float tmp, m, tx, ty;
	register bool flip=0;
	//printf("%i %i\n", std::max(1, diagonal+1-bandwidth), std::min(L_x, diagonal+bandwidth+L_y));
	for( i=std::max(1, diagonal+1-bandwidth); i<std::min(L_x, diagonal+bandwidth+L_y); ++i){
		MM_last    = M[flip];
		II_x_last  = I_x[flip];
		II_y_last  = I_y[flip];
		flip       = !flip;
		MM         = M[flip];
		II_x       = I_x[flip];
		II_y       = I_y[flip];
		MTRX       = matrix[x[i-1]];
		j_offset   = (i-1)-diagonal-bandwidth;
		j_start    = std::max(diagonal-((int)i-1)+bandwidth+1,1);
		j_stop     = std::min(2*bandwidth+1, bandwidth-(int)i+diagonal+(int)L_y);
		for( j=j_start; j<=j_stop; ++j){
			jj   = j+j_offset;
			diag =  i-jj;	
			tmp  = MTRX[y[jj-1]];
			m    = MM_last[j]   + tmp;
			tx   = II_x_last[j] + tmp;
			ty   = II_y_last[j] + tmp;			

			if( (m>=tx) && (m>=ty) ){	
				MM[j]   = m;
				T[i][j] = 1;
			}else if( (tx>=m) && (tx>=ty) ){
				MM[j]   = tx;
				T[i][j] = 0;
			}else{
				MM[j]   = ty;
				T[i][j] = 2;
			}

			if( diag!=lower_diagonal ){
				m  = MM_last[j+1]   - G;
				tx = II_x_last[j+1] - E;	
				if( m>=tx ){
					II_x[j] = m;
					
				}else{
					II_x[j] = tx;
					T[i][j] += 4;
				}
			}else{
				II_x[j] = -INFINITY;
			}

			if( diag!=upper_diagonal ){
				m  = MM[j-1]   - G;
				ty = II_y[j-1] - E;		
				if( m>=ty ){
					II_y[j] = m;
					
				}else{
					II_y[j] = ty;
					T[i][j] += 8;
				}
			}else{
				II_y[j] = -INFINITY;
			}

			if( (jj==L_y-1) || (i==L_x-1) ){
				if( MM[j]>score ){
					max_i = i;
					max_j = jj;
					max_j_band = j;
					score = MM[j];
				}
			}
		}
	}
	//printf("maxi:%i maxj:%i jband:%i score:%f\n", max_i, max_j, max_j_band, score);
	if(verbose) std::cout << "done." << std::endl;
}



void SWaln::_align_local_band_flip(const int diagonal, const int bandwidth) throw (std::exception){
	if(verbose) std::cout << "Fill matrices (LOCAL_BAND) ..." << std::endl;

	if( (diagonal<-((int)L_y-2)) || (diagonal>((int)L_x-2)) ){
		throw MyException("Diagonal(%i) is out of range! Possible values are [%i..%i].\n",diagonal ,-(L_y-2), (L_x-2));
	}
	
	score = -INFINITY;

	//main diagonal is 0   d=i-j!
	const int upper_diagonal = diagonal + bandwidth;
	const int lower_diagonal = diagonal - bandwidth;

	register const float *MTRX;
	register float *MM, *MM_last;
	register float *II_x, *II_x_last;
	register float *II_y, *II_y_last;
	// matrix positions, -- i-1 position int the sequence
	// j position in the band
	// j_start start of band idx
	// j_stop last position of band
	register int i, j, jj, j_start, j_stop, j_offset, diag;
	//tmp holds different temporary values e.g. BLOSUM score of amino acids
	//m: BLOSUM score + last match score
	//tx: BLOSUM score + last score of a gap in x 
	//ty: BLOSUM score + last score of a gap in y
	register float tmp, m, tx, ty;
	register bool flip=0;
	//printf("%i %i\n", std::max(1, diagonal+1-bandwidth), std::min(L_x, diagonal+bandwidth+L_y));
	for( i=std::max(1, diagonal+1-bandwidth); i<std::min(L_x, diagonal+bandwidth+L_y); ++i){
		MM_last    = M[flip];
		II_x_last  = I_x[flip];
		II_y_last  = I_y[flip];
		flip       = !flip;
		MM         = M[flip];
		II_x       = I_x[flip];
		II_y       = I_y[flip];
		MTRX       = matrix[x[i-1]];
		j_offset   = (i-1)-diagonal-bandwidth;
		j_start    = std::max(diagonal-((int)i-1)+bandwidth+1,1);
		j_stop     = std::min(2*bandwidth+1, bandwidth-(int)i+diagonal+(int)L_y);
		for( j=j_start; j<=j_stop; ++j){
			jj   = j+j_offset;
			diag =  i-jj;
			tmp  = MTRX[y[jj-1]];
			m    = MM_last[j]   + tmp;
			tx   = II_x_last[j] + tmp;
			ty   = II_y_last[j] + tmp;			

			if( (m>0) && (m>=tx) && (m>=ty) ){	
				MM[j]   = m;
				T[i][j] = 1;
			}else if( (tx>0) && (tx>=m) && (tx>=ty) ){
				MM[j]   = tx;
				T[i][j] = 0;
			}else if( (ty>0) && (ty>=m) && (ty>=tx) ){
				MM[j]   = ty;
				T[i][j] = 2;
			}else{
				MM[j]   = 0;
				T[i][j] = 3;
			}

			if( diag!=lower_diagonal ){
				m  = MM_last[j+1]   - G;
				tx = II_x_last[j+1] - E;	
				if( m>=tx ){
					II_x[j] = m;
					
				}else{
					II_x[j] = tx;
					T[i][j] += 4;
				}
			}else{
				II_x[j] = -INFINITY;
			}

			if( diag!=upper_diagonal ){
				m  = MM[j-1]   - G;
				ty = II_y[j-1] - E;		
				if( m>=ty ){
					II_y[j] = m;
					
				}else{
					II_y[j] = ty;
					T[i][j] += 8;
				}
			}else{
				II_y[j] = -INFINITY;
			}

			if( MM[j]>score ){
				max_i = i;
				max_j = jj;
				max_j_band = j;
				score = MM[j];
			}
		}
	}
	//printf("maxi:%i maxj:%i jband:%i score:%f\n", max_i, max_j, max_j_band, score);
	if(verbose) std::cout << "done." << std::endl;
}


std::ostream& SWaln::print(std::ostream &out, const size_t width){
	Alignment *a = new Alignment();
	traceback( *a );	
	char buffer[width+1];
//	char mbuffer[width+1];
//	char format[128];
	char xbuffer[width+1];
	char ybuffer[width+1];
	size_t jj, i, j, len_aln;
	len_aln = a->f.len_aln;
	for(i=0; i<len_aln; i=std::min(i+width, len_aln) ){
		jj=0;

		for(j=i; j<std::min(i+width, len_aln); ++j){
			if( a->x[j]==a->y[j] )
				buffer[jj] = '|';
			else if( (a->x[j]<AMINOACID_DIM) && (a->y[j]<AMINOACID_DIM)){
				if( (matrix[a->x[j]][a->y[j]] > 0) )   buffer[jj] = ':';
				else                                   buffer[jj] = '~';
			}else buffer[jj] = ' ';
			if( a->x[j]<AMINOACID_DIM ) xbuffer[jj] = int2aa[ a->x[j] ];
			else xbuffer[jj] = a->x[j];
			if( a->y[j]<AMINOACID_DIM ) ybuffer[jj] = int2aa[ a->y[j] ];
			else ybuffer[jj] = a->y[j];
			++jj;
		}
		//sprintf(format,"\%10
		//sprintf(mbuffer,"%10i%" %10i", i, std::min(i+width, len_aln));
		out.write(xbuffer, jj);	
		out << std::endl;
		out.write(buffer, jj);
		out << std::endl;
		out.write(ybuffer, jj);
		out << std::endl << std::endl << std::endl;
	}
	
	out << "Length of sequence one             : " << a->f.len_x     << std::endl;
	out << "Length of sequence two             : " << a->f.len_y     << std::endl;
	out << "Length of  alignment               : " << a->f.len_aln   << std::endl;
	out << "Equivalenced residues (incl. gaps) : " << a->f.len_equ   << std::endl;
	out << "Score (for equiv. residues)        : " << a->f.score     << std::endl;
	out << "Identities                         : " << a->f.idents    << std::endl;
	out << "Positive matches (incl. idents)    : " << a->f.pos_matches << std::endl;
	out << "Less or equal to zero matches      : " << a->f.neg_zero_matches << std::endl;
	out << "Gaps                               : " << a->f.gaps      << std::endl;
	out << "Percentage identity                : " <<  a->f.idents/(double)(a->f.len_equ) << std::endl;

	if(verbose) _print_memory_usage(out);

	//std::cerr << "REAL SCORE:" << _compute_overlap_score(a->x, a->y, a->f.len_aln) <<std::endl;

	delete [] a->x;
	delete [] a->y;	
	delete a;
	return out;
}


void SWaln::traceback( Alignment &aln ){
	switch(method){
		case SWaln::OVERLAP_BAND:{
			_traceback_overlap_band(aln);
			break;
		}
		case SWaln::LOCAL_FLIP:{
			_traceback_local(aln);
			break;
		}
		case SWaln::LOCAL_BAND:{
			_traceback_local_band(aln);
			break;
		}
		default:{
			_traceback_overlap(aln);
		}
	}
}

Pair_aln* SWaln::traceback_ij(){
	SWaln::Alignment aln;
	switch(method){
		case SWaln::OVERLAP_BAND:{
			_traceback_overlap_band(aln);
			break;
		}
		case SWaln::LOCAL_FLIP:{
			_traceback_local(aln);
			break;
		}
		case SWaln::LOCAL_BAND:{
			_traceback_local_band(aln);
			break;
		}
		default:{
			_traceback_overlap(aln);
		}
	}
	Pair_aln *ret = new Pair_aln(aln.f.len_equ);
	ret->len_x            = aln.f.len_x;
	ret->len_y            = aln.f.len_y;
	ret->len_aln          = aln.f.len_aln;
	ret-> len_equ         = aln.f.len_equ;
	ret->idents           = aln.f.idents;
	ret->gaps             = aln.f.gaps;
	ret->pos_matches      = aln.f.pos_matches;
	ret->neg_zero_matches = aln.f.neg_zero_matches;
	int i=max_i-1;
	int j=max_j-1;
	int alni = aln.f.len_equ-1;
	for( int pos=aln.f.len_equ-1; pos>=0; --pos ){
		if( aln.x[alni]!='-' && aln.y[alni]!='-' ){
			ret->x[pos] = i;
			ret->y[pos] = j;
			--i;--j;
		}else 	if( aln.x[alni]=='-' ){
				ret->x[pos] = -1;
				ret->y[pos] = j;
				--j;
			}else{
				ret->x[pos] = i;
				ret->y[pos] = -1;
				--i;
			}
		--alni;
	}
	delete [] aln.x;
	delete [] aln.y;
	return ret;	
}

void SWaln::_traceback_overlap(Alignment &aln){

	if(verbose) std::cout << "Traceback (OVERLAP)..." << std::endl;
	char *alnx = new char[L_x+L_y];
	char *alny = new char[L_x+L_y];
	size_t aln_pos = L_x+L_y-1;	
	int i,j;
	// print overhanging ends
	i=L_x-1;	
	while(i>max_i){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = '-';
		--aln_pos; --i;
	}

	j=L_y-1;
	while(j>max_j){
		alnx[aln_pos] = '-';
		alny[aln_pos] = y[j-1];
		--aln_pos; --j;
	}

	// traceback!
	size_t equ_matches = aln_pos;
	while( (i>0) && (j>0) ){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = y[j-1];
		--aln_pos;
		--j;--i;
		if( (T[i+1][j+1]&3)==1 ) continue;
		if( (T[i+1][j+1]&3)==0){
			alnx[aln_pos] = x[i-1];
			alny[aln_pos] = '-';
			--aln_pos;

			while( ((T[i][j]&4)>>2)==1 && i>0){
				--i;
				alnx[aln_pos] = x[i-1];
				alny[aln_pos] = '-';
				--aln_pos;
				
			}
			--i;
		}else{ //T[i+1][j+1]==2
			alnx[aln_pos] = '-';
			alny[aln_pos] = y[j-1];
			--aln_pos;
			while( ((T[i][j]&8)>>3)==1 && j>0){
				--j;
				alnx[aln_pos] = '-';
				alny[aln_pos] = y[j-1];
				--aln_pos;
				
			}
			--j;
		}		
	}
	equ_matches = equ_matches - aln_pos;
	
	// print overhanging ends	
	while(i>0){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = '-';
		--aln_pos; --i;
	}
	while(j>0){
		alnx[aln_pos] = '-';
		alny[aln_pos] = y[j-1];
		--aln_pos; --j;
	}

	int aln_len = L_x+L_y-aln_pos-1;
	
	memmove(alnx, &alnx[aln_pos+1], aln_len);
	memmove(alny, &alny[aln_pos+1], aln_len);

	int idents=0, pos_matches=0, neg_matches=0, gaps=0;
	for(i=0; i<aln_len; ++i ){
		if( alnx[i]==alny[i] )
			++idents;
		else if( (alnx[i]<AMINOACID_DIM) && (alny[i]<AMINOACID_DIM))
			if( (matrix[alnx[i]][alny[i]] > 0) )  
				++pos_matches;
			else
				++neg_matches;			
		else ++gaps;
	}

	aln.x             = alnx;
	aln.y             = alny;
	aln.f.gaps        = gaps;
	aln.f.len_equ     = equ_matches;
	aln.f.len_aln     = aln_len;
	aln.f.score       = score;
	aln.f.len_x       = L_x-1;
	aln.f.len_y       = L_y-1;
	aln.f.idents      = idents;
	aln.f.pos_matches = pos_matches+idents;
	aln.f.neg_zero_matches = neg_matches;

	if(verbose) std::cout << "done." << std::endl << std::endl;
}


void SWaln::_traceback_local(Alignment &aln){

	if(verbose) std::cout << "Traceback (LOCAL)..." << std::endl;
	char *alnx = new char[L_x+L_y];
	char *alny = new char[L_x+L_y];
	size_t aln_pos = L_x+L_y-1;	
	int i = max_i;
	int j = max_j;

	// traceback!
	size_t equ_matches = aln_pos;
	while( (i>0) && (j>0) ){
		if( (T[i][j]&3)==3 ) break;
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = y[j-1];
		--aln_pos;
		--j;--i;
		if( (T[i+1][j+1]&3)==1 ) continue;
		if( (T[i+1][j+1]&3)==0){
			alnx[aln_pos] = x[i-1];
			alny[aln_pos] = '-';
			--aln_pos;

			while( ((T[i][j]&4)>>2)==1 && i>0){
				--i;
				alnx[aln_pos] = x[i-1];
				alny[aln_pos] = '-';
				--aln_pos;			
			}
			--i;
		}else{ //T[i+1][j+1]==2
			alnx[aln_pos] = '-';
			alny[aln_pos] = y[j-1];
			--aln_pos;
			while( ((T[i][j]&8)>>3)==1 && j>0){
				--j;
				alnx[aln_pos] = '-';
				alny[aln_pos] = y[j-1];
				--aln_pos;	
			}
			--j;
		}		
	}
	equ_matches = equ_matches - aln_pos;
	
	_xstart = i;
	_ystart = j;

	int aln_len = L_x+L_y-aln_pos-1;
	
	memmove(alnx, &alnx[aln_pos+1], aln_len);
	memmove(alny, &alny[aln_pos+1], aln_len);

	int idents=0, pos_matches=0, neg_matches=0, gaps=0;
	for(i=0; i<aln_len; ++i ){
		if( alnx[i]==alny[i] )
			++idents;
		else if( (alnx[i]<AMINOACID_DIM) && (alny[i]<AMINOACID_DIM))
			if( (matrix[alnx[i]][alny[i]] > 0) )  
				++pos_matches;
			else
				++neg_matches;			
		else ++gaps;
	}

	aln.x             = alnx;
	aln.y             = alny;
	aln.f.gaps        = gaps;
	aln.f.len_equ     = equ_matches;
	aln.f.len_aln     = aln_len;
	aln.f.score       = score;
	aln.f.len_x       = L_x-1;
	aln.f.len_y       = L_y-1;
	aln.f.idents      = idents;
	aln.f.pos_matches = pos_matches+idents;
	aln.f.neg_zero_matches = neg_matches;

	if(verbose) std::cout << "done." << std::endl << std::endl;
}


void SWaln::_traceback_overlap_band(Alignment &aln){
	if(verbose) std::cout << "Traceback (OVERLAP_BAND)..." << std::endl;
	char *alnx = new char[L_x+L_y];
	char *alny = new char[L_x+L_y];
	size_t aln_pos = L_x+L_y-1;	
	int i,j;
	// print overhanging ends
	i=L_x-1;	
	while(i>max_i){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = '-';
		--aln_pos; --i;
	}

	j=L_y-1;
	while(j>max_j){
		alnx[aln_pos] = '-';
		alny[aln_pos] = y[j-1];
		--aln_pos; --j;
	}

	// traceback!
	size_t equ_matches = aln_pos;
	while( (i>0) && (j>0) ){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = y[j-1];
		--aln_pos;
		--i;--j;
		if( (T[i+1][max_j_band]&3)==1 ) continue;
		if( (T[i+1][max_j_band]&3)==0){
			alnx[aln_pos] = x[i-1];
			alny[aln_pos] = '-';
			--aln_pos; 
			while( ((T[i][max_j_band]&4)>>2)==1 && i>0){
				--i;
				alnx[aln_pos] = x[i-1];
				alny[aln_pos] = '-';
				--aln_pos;
				++max_j_band;
			}
			++max_j_band;
			--i;
		}else{ //T[i+1][j+1]==2
			alnx[aln_pos] = '-';
			alny[aln_pos] = y[j-1];
			--aln_pos;
			while( ((T[i][max_j_band]&8)>>3)==1 && j>0){
				--j;
				alnx[aln_pos] = '-';
				alny[aln_pos] = y[j-1];
				--aln_pos;
				--max_j_band;
			}
			--max_j_band;
			--j;
		}		
	}
	equ_matches = equ_matches - aln_pos;
	
	// print overhanging ends	
	while(i>0){
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = '-';
		--aln_pos; --i;
	}
	while(j>0){
		alnx[aln_pos] = '-';
		alny[aln_pos] = y[j-1];
		--aln_pos; --j;
	}

	int aln_len = L_x+L_y-aln_pos-1;
	
	memmove(alnx, &alnx[aln_pos+1], aln_len);
	memmove(alny, &alny[aln_pos+1], aln_len);

	int idents=0, pos_matches=0, neg_matches=0, gaps=0;
	for(i=0; i<aln_len; ++i ){
		if( alnx[i]==alny[i] )
			++idents;
		else if( (alnx[i]<AMINOACID_DIM) && (alny[i]<AMINOACID_DIM))
			if( (matrix[alnx[i]][alny[i]] > 0) )  
				++pos_matches;
			else
				++neg_matches;			
		else ++gaps;
	}

	aln.x             = alnx;
	aln.y             = alny;
	aln.f.gaps        = gaps;
	aln.f.len_equ     = equ_matches;
	aln.f.len_aln     = aln_len;
	aln.f.score       = score;
	aln.f.len_x       = L_x-1;
	aln.f.len_y       = L_y-1;
	aln.f.idents      = idents;
	aln.f.pos_matches = pos_matches+idents;
	aln.f.neg_zero_matches = neg_matches;

	if(verbose) std::cout << "done." << std::endl << std::endl;
}


void SWaln::_traceback_local_band(Alignment &aln){
	if(verbose) std::cout << "Traceback (LOCAL_BAND)..." << std::endl;
	char *alnx = new char[L_x+L_y];
	char *alny = new char[L_x+L_y];
	size_t aln_pos = L_x+L_y-1;	
	int i = max_i;
	int j = max_j;
	//std::cerr << "i:" << i << " j:" << j << std::endl;
	// traceback!
	size_t equ_matches = aln_pos;
	while( (i>0) && (j>0) ){
		if( (T[i][max_j_band]&3)==3 ) break;
		alnx[aln_pos] = x[i-1];
		alny[aln_pos] = y[j-1];
		--aln_pos;
		--i;--j;
		if( (T[i+1][max_j_band]&3)==1 ) continue;
		if( (T[i+1][max_j_band]&3)==0){
			alnx[aln_pos] = x[i-1];
			alny[aln_pos] = '-';
			--aln_pos; 
			while( ((T[i][max_j_band]&4)>>2)==1 && i>0){
				--i;
				alnx[aln_pos] = x[i-1];
				alny[aln_pos] = '-';
				--aln_pos;
				++max_j_band;
			}
			++max_j_band;
			--i;
		}else{ //T[i+1][j+1]==2
			alnx[aln_pos] = '-';
			alny[aln_pos] = y[j-1];
			--aln_pos;
			while( ((T[i][max_j_band]&8)>>3)==1 && j>0){
				--j;
				alnx[aln_pos] = '-';
				alny[aln_pos] = y[j-1];
				--aln_pos;
				--max_j_band;
			}
			--max_j_band;
			--j;
		}		
	}
	equ_matches = equ_matches - aln_pos;
	
	int aln_len = L_x+L_y-aln_pos-1;
	
	memmove(alnx, &alnx[aln_pos+1], aln_len);
	memmove(alny, &alny[aln_pos+1], aln_len);

	int idents=0, pos_matches=0, neg_matches=0, gaps=0;
	for(i=0; i<aln_len; ++i ){
		if( alnx[i]==alny[i] )
			++idents;
		else if( (alnx[i]<AMINOACID_DIM) && (alny[i]<AMINOACID_DIM))
			if( (matrix[alnx[i]][alny[i]] > 0) )  
				++pos_matches;
			else
				++neg_matches;			
		else ++gaps;
	}

	aln.x             = alnx;
	aln.y             = alny;
	aln.f.gaps        = gaps;
	aln.f.len_equ     = equ_matches;
	aln.f.len_aln     = aln_len;
	aln.f.score       = score;
	aln.f.len_x       = L_x-1;
	aln.f.len_y       = L_y-1;
	aln.f.idents      = idents;
	aln.f.pos_matches = pos_matches+idents;
	aln.f.neg_zero_matches = neg_matches;

	if(verbose) std::cout << "done." << std::endl << std::endl;
}


void SWaln::_init_overlap_matrices_full() throw (std::exception){
	if(verbose) std::cout << "Allocating memory (OVERLAP_FULL) ..." << std::endl;
	M   = new  float*[L_x];
	I_x = new  float*[L_x];
	I_y = new  float*[L_x];
	T   = new char*[L_x];

	for( int i=0; i<L_x; ++i ){
		M[i]   = new  float[L_y];
		I_x[i] = new  float[L_y];
		I_y[i] = new  float[L_y];
		T[i]   = new char[L_y];
	}

	//init matrices
	for( int i=1; i<L_x; ++i ){
			I_y[i][0] = 0 ;//-(GAP_OPEN + (i*GAP_EXTENSION));
			I_x[i][0] = -INFINITY;
			  M[i][0] = 0;// -(GAP_OPEN + (i*GAP_EXTENSION));
			  T[i][0] = 0;
	}
	for( int i=1; i<L_y; ++i ){
			I_x[0][i] = 0 ;//-(GAP_OPEN + (j*GAP_EXTENSION));
			I_y[0][i] = -INFINITY;
			  M[0][i] = 0;//-(GAP_OPEN + (j*GAP_EXTENSION));
			  T[0][i] = 0;
  	}	

	I_x[0][0] = -INFINITY;
	I_y[0][0] = -INFINITY;
	M[0][0]   = 0;
	T[0][0]   = 0;

	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_delete_overlap_matrices_full(){
	if(verbose) std::cout << "Deleting allocated memory (OVERLAP_FULL) ..." << std::endl;
	for( int i=0; i<L_x; ++i ){
		delete [] M[i];
		delete [] I_x[i];
		delete [] I_y[i];
		delete [] T[i];
	}

	delete [] M;
	delete [] I_x;
	delete [] I_y;
	delete [] T;
	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_init_overlap_matrices_flip() throw (std::exception){
	if(verbose) std::cout << "Allocating memory (OVERLAP_FLIP) ..." << std::endl;
	M   = new  float*[2];
	I_x = new  float*[2];
	I_y = new  float*[2];
	T   = new char*[L_x];


	for(int i=0; i<L_x; ++i ){
		if( i<2){
			M[i]   = new  float[L_y];
			I_x[i] = new  float[L_y];
			I_y[i] = new  float[L_y];
		}
		T[i]   = new char[L_y];
	}

	//init matrices
	for(int i=1; i<L_x; ++i ){
			if(i<2){
				I_y[i][0] = -INFINITY;
				I_x[i][0] = -INFINITY;
			 	 M[i][0] = 0;
			}
			T[i][0] = 0;
	}

	for(int i=1; i<L_y; ++i ){
			I_x[0][i] = -INFINITY;//-(GAP_OPEN + (j*GAP_EXTENSION));
			I_y[0][i] = -INFINITY;
			  M[0][i] = 0;//-(GAP_OPEN + (j*GAP_EXTENSION));
			  T[0][i] = 0;
  	}	

	I_x[0][0] = -INFINITY;
	I_y[0][0] = -INFINITY;
	M[0][0]   = 0;
	T[0][0]   = 0;

	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_delete_overlap_matrices_flip(){
	if(verbose) std::cout << "Deleting allocated memory (OVERLAP_FLIP) ..." << std::endl;
	for(int i=0; i<L_x; ++i ){
		if(i<2){
			delete [] M[i];
			delete [] I_x[i];
			delete [] I_y[i];
		}
		delete [] T[i];
	}

	delete [] M;
	delete [] I_x;
	delete [] I_y;
	delete [] T;
	if(verbose) std::cout << "done." << std::endl;
}

void SWaln::_init_overlap_matrices_band_flip() throw (std::exception){
	if(verbose) std::cout << "Allocating memory (OVERLAP_BAND, LOCAL_BAND) ..." << std::endl;
	M   = new  float*[2];
	I_x = new  float*[2];
	I_y = new  float*[2];
	T   = new char*[L_x];

	int total_bandwidth = 2 + (2*W); 
	M[0]   = new  float[total_bandwidth];
	I_x[0] = new  float[total_bandwidth];
	I_y[0] = new  float[total_bandwidth];
	M[1]   = new  float[total_bandwidth];
	I_x[1] = new  float[total_bandwidth];
	I_y[1] = new  float[total_bandwidth];


	for(int i=0; i<L_x; ++i )	T[i] = new char[total_bandwidth];

	//init matrices
	for(int j=0; j<total_bandwidth; ++j ){
		I_x[0][j] = -INFINITY;
		I_y[0][j] = -INFINITY;
	  	  M[0][j] = 0;
		I_x[1][j] = -INFINITY;
		I_y[1][j] = -INFINITY;
	  	  M[1][j] = 0;
	}

	I_x[0][0] = -INFINITY;
	I_y[0][0] = -INFINITY;
	M[0][0] = 0;
	T[0][0] = 0;

	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::_delete_overlap_matrices_band_flip(){	
	if(verbose) std::cout << "Deleting allocated memory (OVERLAP_BAND, LOCAL_BAND) ..." << std::endl;

	for(int i=0; i<L_x; ++i ){
		if( i<2 ){
			delete [] M[i];
			delete [] I_x[i];
			delete [] I_y[i];
		}
		delete [] T[i];
	}
	delete [] M;
	delete [] I_x;
	delete [] I_y;
	delete [] T;
	if(verbose) std::cout << "done." << std::endl;
}


void SWaln::compute_features( Features& f){
	if(verbose) std::cout << "Computing features ..." << std::endl;
	f.len_x            = L_x-1;
	f.len_y            = L_y-1;
	f.idents           = 0;
	f.len_aln          = 0;
	f.len_equ          = 0;
	f.score            = 0;
	f.gaps             = 0;
	f.neg_zero_matches = 0;
	f.pos_matches      = 0;	

	size_t aln_pos = L_x+L_y-1;	
	
	f.gaps  += L_x-1-max_i;
	f.gaps  += L_y-1-max_j;
	aln_pos -= L_x-1-max_i;
	aln_pos -= L_y-1-max_j;

	int i=max_i, j=max_j;
	f.len_equ = aln_pos;
	while( (i>0) && (j>0) ){
		if( x[i-1]==y[j-1] ) ++f.idents;
		if( matrix[x[i-1]][y[j-1]]>0 ) ++f.pos_matches;
		else ++f.neg_zero_matches;
		--aln_pos;
		--j;--i;
		if( (T[i+1][j+1]&3)==1 ) continue;
		if( (T[i+1][j+1]&3)==0){
			--aln_pos;++f.gaps;
			while( ((T[i][j]&4)>>2)==1 ){
				--i;
				--aln_pos;
				++f.gaps;			
			}
			--i;
		}else{ //T[i+1][j+1]==2
			--aln_pos;++f.gaps;
			while( ((T[i][j]&8)>>3)==1 ){
				--j;
				--aln_pos;
				++f.gaps;			
			}
			--j;
		}		
	}

	f.gaps  += i+j;
	f.len_equ = f.len_equ - aln_pos;
	aln_pos   = aln_pos-i-j;	
	f.len_aln = L_x+L_y-aln_pos-1;
	f.score   = score;
	f.len_x   = L_x-1;
	f.len_y   = L_y-1;
	
	if(verbose) std::cout << "done." << std::endl << std::endl;
}


std::ostream& SWaln::_print_memory_usage(std::ostream& out){
	std::string fn = "/proc/self/status";
	std::ifstream in( fn.c_str() );
	if(!in){
		out << "Cannot open " << fn << std::endl;
		return out;
	}
	
	out << std::endl << "MEMORY:" << std::endl;
	std::string key;
	while(getline(in,key)){
		if( key.find("VmPeak")!=std::string::npos )
			out << key << std::endl;
		else if( key.find("VmSize")!=std::string::npos )
			out << key << std::endl;
		else if( key.find("VmHWM")!=std::string::npos )
			out << key << std::endl;
		else if( key.find("VmRSS")!=std::string::npos )
			out << key << std::endl;
		else if( key.find("VmData")!=std::string::npos )
			out << key << std::endl;

	}
	in.close();
	return out;
}

void SWaln::_print_matrix(char **MM){
	printf("\n\nM:\n");

	for(int j=0; j<L_y; ++j){
		if(j==0){
			printf("%6c ", ' ');
			printf("%6c ", ' ');
			for(int i=1; i<L_x; ++i){
				printf("%8c    ", int2aa[x[i-1]]);
			}
			printf("\n");
		}
		for(int i=0; i<L_x; ++i){			
			if( (i==0) && (j>0) ) printf("%8c    ", int2aa[y[j-1]]);
			else if(i==0) printf("%6c ", ' ');
			printf("%2i,%2i(%i,%i,%i) ", i,j,T[i][j]&3,T[i][j]&4 >> 2,T[i][j]&8>>3);
			//printf("%+3i(%i,%i,%i) ", MM[i][j], T[i][j]&3,T[i][j]&4 >> 2,T[i][j]&8>>3);
		}
		printf("\n");
	}
}



float SWaln::_compute_overlap_score(char *alnx, char *alny, size_t len){
	
	float ret=0;
	
	size_t start=0;
	for(size_t i=0; i<len; ++i){
		if( (alnx[i]=='-') || (alny[i]=='-') ) ++start;
		else break;
	}
	size_t stop=len;
	for(size_t i=len-1; i>=0; --i){
		if( (alnx[i]=='-') || (alny[i]=='-') ) --stop;
		else break;
	}
	//std::cout << start << " " <<stop << std::endl;
	bool open = false;
	for(size_t i=start; i<stop; ++i){
		//std::cout << "RET:" << ret << std::endl;
		if(alnx[i]=='-' || alny[i]=='-'){
			if(!open) ret -= G;
			else ret -= E;
			open = true;
		}else{ 
			//std::cout << matrix[ alnx[i] ][ alny[i] ] << std::endl;
			ret += matrix[ alnx[i] ][ alny[i] ];
			open=false;
		}
	}

	return ret;
}

void SWaln::print_traceback( Alignment &aln ){
	size_t l=_xstart, k=_ystart;
	for( int i=0; i<aln.f.len_aln; ++i ){
		std::cout << l << " " << k << std::endl;
		if( aln.x[i]!='-' ) ++l;
		if( aln.y[i]!='-' ) ++k;
	} 
}

