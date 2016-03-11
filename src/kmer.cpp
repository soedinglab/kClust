/***************************************************************************
 *   Copyright (C) 2006 by christian mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "kmer.h"

Kmer::Kmer(const size_t k, const double thr, bool score_correction, Kmer_matrix* m_4mer, Kmer_matrix* m_2mer):
		k(k), 
		thr(thr),
		k_thr(k*thr),
		sc(score_correction),
		m_4mer(m_4mer),
		m_2mer(m_2mer){
	
}

void Kmer::print_scoring_matrix(size_t imax, size_t jmax){
	for (int i = 0; i < imax; i++){
		for (int j = 0; j < jmax; j++){
			std::cout << aa_matrix[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void Kmer::init_to_sequence(Sequence * s){
	Matrix * m = s->get_matrix();
	AMINOACID_DIM = m->get_aa_dim();
	int2aa = m->get_int2aa();
	aa2int = m->get_aa2int();
	aa_powers = m->get_aa_powers();
	_dimers = m->get_dimers();
	
	aa_matrix = s->get_scoring_matrix();
	max_scores = s->get_max_scores();
	score_corr = s->get_score_correction();
}


Kmer::~Kmer(){
}

size_t Kmer::get_memory_usage(){
	size_t ret = 3*k;
	ret += seed_list.get_instance_memory_usage();
	return ret;	
}

Kmer_seed_list::Iterator Kmer::begin() { return seed_list.begin();	   }
Kmer_seed_list::Iterator Kmer::end()   { return seed_list.end();	      }
size_t Kmer::get_list_length()         { return seed_list.load();       }
size_t Kmer::get_kmer_count()          { return kmer_count;             }

inline void Kmer::_init_smax(int* kmer, float* smax, float* corr){
	float sum=0;
	for(int i=k-1; i>=0; --i){
		smax[i] = k_thr-sum;
		sum += max_scores[kmer[i]] + corr[i];
	}
}

inline void Kmer::_naiv_kmer_sort(int* kmer, int* sort){
	int tmp=0;
	int *is=kmer;
	int *ptr =0;
	int *ie = &kmer[k];	
	int *sis;
	int *sptr;
	for( size_t i=0; i<k; ++i ) sort[i]=i;	
	sis = sort;	
	while( is!=ie ){
		ptr  = is;
		sptr = sis;
		while( ptr!=ie ){
			if( max_scores[*ptr]>max_scores[*is] ){
				tmp     = *is;
				*is     = *ptr;
				*ptr    = tmp;
				tmp     = *sis;
				*sis    = *sptr;
				*sptr   = tmp;
			}
			++ptr;
			++sptr;
		}
		++is;
		++sis;
	}
}

void Kmer::_sort(float* corr, int* sort){
	float tmp[k];
	memcpy(&tmp, corr, k*sizeof(float));
	for (int i = 0; i < k; i++) corr[i] = tmp[sort[i]];
}


size_t Kmer::create_kmer_list(Sequence *q){
	int kmer[k];
	int sort[k];
	float smax[k];
	float corr[k];
	
	init_to_sequence(q);
	seed_list.reset();
	if( q->length()<k ) return 0;
	const size_t Lk         = q->length()-k;
	kmer_count              = 0;	
	for( size_t i=0; i<=Lk; ++i ){
		q->get_kmer_at(kmer, k, i);
		memcpy(&corr, &score_corr[i], k*sizeof(float));
		_naiv_kmer_sort(kmer, sort);
		_sort(corr, sort);
		_init_smax(kmer, smax, corr);
		_add_similar_kmers(kmer, smax, sort, i, corr);
	}
	return kmer_count;	
}

size_t Kmer::create_kmer_list_spaced(Sequence *q) throw (std::exception){

	int kmer[k];
	int sort[k];
	float smax[k];
	float corr[k];
	memset(corr, 0.0f, k*sizeof(float));
	
	init_to_sequence(q);
	seed_list.reset();
	if( k==6 && q->length()<(k+2) ) return 0;
	if( q->length()<(k+1) )         return 0;
	if(k<3) throw MyException("Spaced k-mer generation implemented for k>2 only!");

	size_t Lk   = q->length()-k-1;
	if(k==6) Lk = q->length()-k-2;

	kmer_count              = 0;	
	for( size_t i=0; i<=Lk; ++i ){
		if(k==6){
			kmer[0] = q->get_kmer_entry_at(i);
			kmer[1] = q->get_kmer_entry_at(i+1);
			kmer[2] = q->get_kmer_entry_at(i+3);
			kmer[3] = q->get_kmer_entry_at(i+4);
			kmer[4] = q->get_kmer_entry_at(i+5);
			kmer[5] = q->get_kmer_entry_at(i+7);
		}else{
			kmer[0] = q->get_kmer_entry_at(i);
			kmer[1] = q->get_kmer_entry_at(i+1);
			kmer[2] = q->get_kmer_entry_at(i+3);
			kmer[3] = q->get_kmer_entry_at(i+4);
		}
		_naiv_kmer_sort(kmer, sort);
		_init_smax(kmer, smax, corr);
		_add_similar_kmers(kmer, smax, sort, i, corr);
	}
	return kmer_count;	
}


void Kmer::_add_similar_kmers(int* kmer, float* smax, int* sort, size_t i, float* corr){
	int kmer0, kmer1, kmer2, kmer3, kmer4, kmer5, kmer6;
	float smax0,  smax1, smax2, smax3, smax4, smax5, smax6;
	switch(k){
		case 3:{
			kmer0    = kmer[0];
			kmer1    = kmer[1];
			kmer2    = kmer[2];
			smax0    = smax[0];
			smax1    = smax[1];
			smax2    = smax[2];
			_add_similar_3mers(i, kmer0, kmer1, kmer2, smax0, smax1, smax2, sort, corr);
			break;
		}
		case 4:{
			kmer0    = kmer[0];
			kmer1    = kmer[1];
			kmer2    = kmer[2];
			kmer3    = kmer[3];
			smax0    = smax[0];
			smax1    = smax[1];
			smax2    = smax[2];
			smax3    = smax[3];
			_add_similar_4mers(i, kmer0, kmer1, kmer2, kmer3, smax0, smax1, smax2, smax3, sort, corr);
			break;
		}
		case 5:{
			kmer0    = kmer[0];
			kmer1    = kmer[1];
			kmer2    = kmer[2];
			kmer3    = kmer[3];
			kmer4    = kmer[4];
			smax0    = smax[0];
			smax1    = smax[1];
			smax2    = smax[2];
			smax3    = smax[3];
			smax4    = smax[4];
			_add_similar_5mers(i, kmer0, kmer1, kmer2, kmer3, kmer4, smax0, smax1, smax2, smax3, smax4, sort, corr);
			break;
		}	
		case 6:{
			kmer0    = kmer[0];
			kmer1    = kmer[1];
			kmer2    = kmer[2];
			kmer3    = kmer[3];
			kmer4    = kmer[4];
			kmer5    = kmer[5];
			smax0    = smax[0];
			smax1    = smax[1];
			smax2    = smax[2];
			smax3    = smax[3];
			smax4    = smax[4];
			smax5    = smax[5];
			_add_similar_6mers(i, kmer0, kmer1, kmer2, kmer3, kmer4, kmer5, smax0, smax1, smax2, smax3, smax4, smax5, sort, corr);
			break;
		}
		case 7:{
			kmer0    = kmer[0];
			kmer1    = kmer[1];
			kmer2    = kmer[2];
			kmer3    = kmer[3];
			kmer4    = kmer[4];
			kmer5    = kmer[5];
			kmer6    = kmer[6];
			smax0    = smax[0];
			smax1    = smax[1];
			smax2    = smax[2];
			smax3    = smax[3];
			smax4    = smax[4];
			smax5    = smax[5];
			smax6    = smax[6];
			_add_similar_7mers(i, kmer0, kmer1, kmer2, kmer3, kmer4, kmer5, kmer6, smax0, smax1, smax2, smax3, smax4, smax5, smax6, sort, corr);
			break;
		}
		default:{
			_rec_add_similar_kmers(0,0,0,i, kmer, sort, smax);
			break;
		}
	}
}

/* j:       current pos in the generated kmer
 * s_j:     score up to position j-1
 * idx:     numerical code of the generated kmer so far
 * pos:     position of the kmer in the sequence
 */
void Kmer::_rec_add_similar_kmers( const size_t j, const float s_j, size_t idx, const size_t pos, int* kmer, int* sort, float* smax){
	register float tmp_score = 0;
	register size_t tmp_idx  = 0;
	for( size_t a=0; a<AMINOACID_DIM; ++a ){
		tmp_score = s_j + aa_matrix[kmer[j]][a];
		if( tmp_score < smax[j] ) continue;
		tmp_idx = idx + (a*aa_powers[sort[j]]);
		if( j<(k-1) ){ 
			_rec_add_similar_kmers(j+1, tmp_score, tmp_idx, pos, kmer, sort, smax);
		}else{
			++kmer_count;
			seed_list.append(tmp_idx, tmp_score, pos);
		}	
	}
}

void Kmer::_add_similar_3mers(const size_t i, int kmer0, int kmer1, int kmer2, float smax0, float smax1, float smax2, int* sort, float* corr){
	register size_t a0,a1,a2;
	size_t index0, index1, index2;
	size_t pow0, pow1, pow2;
	register float score0, score1, score2;
	register const float *AA0, *AA1, *AA2;
	
	pow0 = aa_powers[sort[0]];
	pow1 = aa_powers[sort[1]];
	pow2 = aa_powers[sort[2]];

	AA0 = aa_matrix[kmer0];

	for( a0=0; a0<AMINOACID_DIM; ++a0 ){	
		score0 = *(AA0++) + corr[0];
		if( score0 < smax0 ) continue;
		AA1 = aa_matrix[kmer1];
		index0 = a0*pow0;
		for( a1=0; a1<AMINOACID_DIM; ++a1 ){
			score1 = score0 + *(AA1++) + corr[1];
			if( score1 < smax1 ) continue;
			AA2 = aa_matrix[kmer2];
			index1 = index0 + a1*pow1;
			for( a2=0; a2<AMINOACID_DIM; ++a2 ){
				score2 = score1 + *(AA2++) + corr[2];
				if( score2 < smax2 ) continue;
				index2 = index1 + a2*pow2;
				++kmer_count;
				seed_list.append(index2, score2, i);
	}	}	}		
}

void Kmer::_add_similar_4mers(const size_t i, int kmer0, int kmer1, int kmer2, int kmer3, float smax0, float smax1, float smax2, float smax3, int* sort, float* corr){
	register size_t a0,a1,a2,a3;
	size_t index0, index1, index2, index3;
	size_t pow0, pow1, pow2, pow3;
	register float score0, score1, score2, score3;
	register const float *AA0, *AA1, *AA2, *AA3;
	
	pow0 = aa_powers[sort[0]];
	pow1 = aa_powers[sort[1]];
	pow2 = aa_powers[sort[2]];
	pow3 = aa_powers[sort[3]];

	AA0 = aa_matrix[kmer0];

	for( a0=0; a0<AMINOACID_DIM; ++a0 ){	
		score0 = *(AA0++) + corr[0];
		if( score0 < smax0 ) continue;
		AA1 = aa_matrix[kmer1];
		index0 = a0*pow0;
		for( a1=0; a1<AMINOACID_DIM; ++a1 ){
			score1 = score0 + *(AA1++) + corr[1];
			if( score1 < smax1 ) continue;
			AA2 = aa_matrix[kmer2];
			index1 = index0 + a1*pow1;
			for( a2=0; a2<AMINOACID_DIM; ++a2 ){
				score2 = score1 + *(AA2++) + corr[2];
				if( score2 < smax2 ) continue;
				AA3 = aa_matrix[kmer3];
				index2 = index1 + a2*pow2;
				for( a3=0; a3<AMINOACID_DIM; ++a3 ){
					score3 = score2 + *(AA3++) + corr[3];
					if( score3 < smax3 ) continue;
//					AA4 = aa_matrix[kmer4];
					index3 = index2 + a3*pow3;
					++kmer_count;
					seed_list.append(index3, score3, i);
	}	}	}	}	
}


void Kmer::_add_similar_5mers(const size_t i, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, float smax0, float smax1, float smax2, float smax3, float smax4, int* sort, float* corr){
	register size_t a0,a1,a2,a3,a4;
	size_t index0, index1, index2, index3, index4;
	size_t pow0, pow1, pow2, pow3, pow4;
	register float score0, score1, score2, score3, score4;
	register const float *AA0, *AA1, *AA2, *AA3, *AA4;

	pow0 = aa_powers[sort[0]];
	pow1 = aa_powers[sort[1]];
	pow2 = aa_powers[sort[2]];
	pow3 = aa_powers[sort[3]];
	pow4 = aa_powers[sort[4]];

	AA0 = aa_matrix[kmer0];

	for( a0=0; a0<AMINOACID_DIM; ++a0 ){	
		score0 = *(AA0++) + corr[0];
		if( score0 < smax0 ) continue;
		AA1 = aa_matrix[kmer1];
		index0 = a0*pow0;
		for( a1=0; a1<AMINOACID_DIM; ++a1 ){
			score1 = score0 + *(AA1++) + corr[1];
			if( score1 < smax1 ) continue;
			AA2 = aa_matrix[kmer2];
			index1 = index0 + a1*pow1;
			for( a2=0; a2<AMINOACID_DIM; ++a2 ){
				score2 = score1 + *(AA2++) + corr[2];
				if( score2 < smax2 ) continue;
				AA3 = aa_matrix[kmer3];
				index2 = index1 + a2*pow2;
				for( a3=0; a3<AMINOACID_DIM; ++a3 ){
					score3 = score2 + *(AA3++) + corr[3];
					if( score3 < smax3 ) continue;
					AA4 = aa_matrix[kmer4];
					index3 = index2 + a3*pow3;
					for( a4=0; a4<AMINOACID_DIM; ++a4 ){
						score4 = score3 + *(AA4++) + corr[4];
						if( score4 < smax4 ) continue;
						index4 = index3 +a4*pow4;
						++kmer_count;
						seed_list.append(index4, score4, i);
	}	}	}	}	}	
}

void Kmer::_add_similar_6mers(const size_t i, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, float smax0, float smax1, float smax2, float smax3, float smax4, float smax5, int* sort, float* corr){
	register size_t a0,a1,a2,a3,a4,a5;
	size_t index0, index1, index2, index3, index4, index5;
	size_t pow0, pow1, pow2, pow3, pow4, pow5;
	register float score0, score1, score2, score3, score4, score5;
	register const float *AA0, *AA1, *AA2, *AA3, *AA4, *AA5;
	
	pow0 = aa_powers[sort[0]];
	pow1 = aa_powers[sort[1]];
	pow2 = aa_powers[sort[2]];
	pow3 = aa_powers[sort[3]];
	pow4 = aa_powers[sort[4]];
	pow5 = aa_powers[sort[5]];

	//spalte der matrix fuer die position in der sequenz an der ersten position des kmers
	AA0 = aa_matrix[kmer0];

	for( a0=0; a0<AMINOACID_DIM; ++a0 ){	
		score0 = *(AA0++) + corr[0];
		if( score0 < smax0 ) continue;
		AA1 = aa_matrix[kmer1];
		index0 = a0*pow0;
		for( a1=0; a1<AMINOACID_DIM; ++a1 ){
			score1 = score0 + *(AA1++) + corr[1];
			if( score1 < smax1 ) continue;
			AA2 = aa_matrix[kmer2];
			index1 = index0 + a1*pow1;
			for( a2=0; a2<AMINOACID_DIM; ++a2 ){
				score2 = score1 + *(AA2++) + corr[2];
				if( score2 < smax2 ) continue;
				AA3 = aa_matrix[kmer3];
				index2 = index1 + a2*pow2;
				for( a3=0; a3<AMINOACID_DIM; ++a3 ){
					score3 = score2 + *(AA3++) + corr[3];
					if( score3 < smax3 ) continue;
					AA4 = aa_matrix[kmer4];
					index3 = index2 + a3*pow3;
					for( a4=0; a4<AMINOACID_DIM; ++a4 ){
						score4 = score3 + *(AA4++) + corr[4];
						if( score4 < smax4 ) continue;
						AA5 = aa_matrix[kmer5];
						index4 = index3 +a4*pow4;
						for( a5=0; a5<AMINOACID_DIM; ++a5 ){
							score5 = score4 + *(AA5++) + corr[5];
							if( score5 < smax5 ) continue;
							index5 = index4 + a5*pow5;
							++kmer_count;
							seed_list.append(index5, score5, i);
	}	}	}	}	}	}	
}

void Kmer::_add_similar_7mers(const size_t i, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, int kmer6, float smax0, float smax1, float smax2, float smax3, float smax4, float smax5, float smax6, int* sort, float* corr){
	register size_t a0,a1,a2,a3,a4,a5,a6;
	size_t index0, index1, index2, index3, index4, index5, index6;
	size_t pow0, pow1, pow2, pow3, pow4, pow5, pow6;
	register float score0, score1, score2, score3, score4, score5, score6;
	register const float *AA0, *AA1, *AA2, *AA3, *AA4, *AA5, *AA6;

	pow0 = aa_powers[sort[0]];
	pow1 = aa_powers[sort[1]];
	pow2 = aa_powers[sort[2]];
	pow3 = aa_powers[sort[3]];
	pow4 = aa_powers[sort[4]];
	pow5 = aa_powers[sort[5]];
	pow6 = aa_powers[sort[6]];

	AA0 = aa_matrix[kmer0];

	for( a0=0; a0<AMINOACID_DIM; ++a0 ){	
		score0 = *(AA0++) + corr[0];
		if( score0 < smax0 ) continue;
		AA1 = aa_matrix[kmer1];
		index0 = a0*pow0;
		for( a1=0; a1<AMINOACID_DIM; ++a1 ){
			score1 = score0 + *(AA1++) + corr[1];
			if( score1 < smax1 ) continue;
			AA2 = aa_matrix[kmer2];
			index1 = index0 + a1*pow1;
			for( a2=0; a2<AMINOACID_DIM; ++a2 ){
				score2 = score1 + *(AA2++) + corr[2];
				if( score2 < smax2 ) continue;
				AA3 = aa_matrix[kmer3];
				index2 = index1 + a2*pow2;
				for( a3=0; a3<AMINOACID_DIM; ++a3 ){
					score3 = score2 + *(AA3++) + corr[3];
					if( score3 < smax3 ) continue;
					AA4 = aa_matrix[kmer4];
					index3 = index2 + a3*pow3;
					for( a4=0; a4<AMINOACID_DIM; ++a4 ){
						score4 = score3 + *(AA4++) + corr[4];
						if( score4 < smax4 ) continue;
						AA5 = aa_matrix[kmer5];
						index4 = index3 +a4*pow4;
						for( a5=0; a5<AMINOACID_DIM; ++a5 ){
							score5 = score4 + *(AA5++) + corr[5];
							if( score5 < smax5 ) continue;
							AA6 = aa_matrix[kmer6];
							index5 = index4 +a5*pow5;
							for( a6=0; a6<AMINOACID_DIM; ++a6 ){
								score6 = score5 + *(AA6++) + corr[6];
								if( score6 < smax6 ) continue;
								index6 = index5 + a6*pow6;
								++kmer_count;
								seed_list.append(index6, score6, i);
	}	}	}	}	}	}	}	
}

/*void Kmer::print_kmer_seed_list(){
	printf("LOAD:%Zu\n", seed_list.load());
	Kmer_seed_list::Iterator it = seed_list.begin();
	while(it!=seed_list.end()){	
		print_kmer((*it)->index, std::cout);
		printf(" (%-8i, %-3.2f, %-5i)\n", (*it)->index, (*it)->score, (*it)->position);
		++it;
	}
}*/

/*void Kmer::_print_kmer(char *ptr){
	for(size_t i=0; i<k; ++i)printf("%c", int2aa[(int)ptr[i]]);
}

std::ostream& Kmer::print_kmer(size_t index, std::ostream &out){
	for(size_t i=0; i<k; ++i){
		size_t aa = index%21;
		out << int2aa[aa];
		index = index-aa;
		index /= AMINOACID_DIM;
	}
	return out;
}*/

void Kmer::index2kmer(size_t index, char *kmer){
	for(size_t i=0; i<k; ++i){
		size_t aa = index%21;
		kmer[i] = aa;
		index = index-aa;
		index /= AMINOACID_DIM;
	}
}

void Kmer::get_kmer_score_with(Sequence *s, Kmer::kmer_score &ret){
	if( s->length()<k ){
		ret.count=0;
		ret.score=0;
		return;
	}
	const size_t Lk = s->length()-k;
	const char * const ptr       = s->get_sequence();
	//put all kmers of target sequence into hash
	Simple_hash<int> h(Lk);
	for( size_t j=0; j<=Lk; ++j ){ 
		h.put(&ptr[j], k, j);
	}

	//sum up max-kmer-score per position in query
	Kmer_seed_list::Iterator it = seed_list.begin();
	const Kmer_seed_list::Iterator end = seed_list.end();
	ret.count   = 0;
	ret.score  = 0.0f;
	float max_kmer_score = 0.0f;
	int last_pos = 0;
	if( it!=end ) last_pos = (*it)->position; 
	
	char mer[k];
	while( it!=end ){
		if( last_pos != (*it)->position ){
			ret.score += max_kmer_score;
			if(max_kmer_score>0.0f) ++ret.count;
			max_kmer_score = 0.0f;
			last_pos = (*it)->position;	
		}
		index2kmer((*it)->index, mer);
		//print_kmer( (*it)->index, std::cerr);
		if( h.exists(mer, k) ){
			max_kmer_score = std::max(max_kmer_score, (*it)->score);
		}
		++it;
	}
	ret.score += max_kmer_score;
	if(max_kmer_score>0.0f) ++ret.count;
}

std::ostream& Kmer::print_kmer_matches_with(Sequence *s, std::ostream &out){
	if( s->length()<k )return out;
	const size_t Lk = s->length()-k;
	const char * const ptr       = s->get_sequence();
	//put all kmers of template sequence into hash
	Simple_hash<int> h(Lk);
	for( size_t j=0; j<=Lk; ++j )h.put(&ptr[j], k, j);

	//sum up max-kmer-score per position in query
	Kmer_seed_list::Iterator it = seed_list.begin();
	const Kmer_seed_list::Iterator end = seed_list.end();
	float max_kmer_score = 0.0f;
	int max_kmer_pos = 0;
	int last_pos = 0;
	if( it!=end ) last_pos = (*it)->position; 
	
	char mer[k];
	while( it!=end ){
		index2kmer((*it)->index, mer);
		if( h.exists(mer, k) ){
			if( max_kmer_score< (*it)->score ){
				max_kmer_score = (*it)->score;
				max_kmer_pos   = (*it)->position;	
			}
		}
		if( last_pos != (*it)->position ){
			if(max_kmer_score>0.0f){
				out << last_pos << " " << max_kmer_pos << std::endl;
			}
			max_kmer_score = 0.0f;
			last_pos = (*it)->position;	
		}
		++it;
	}
	return out;
}

//void Kmer::get_kmer_score_spaced_with(Sequence *s, Kmer::kmer_score &ret){
	// example PATTERN: XX0XXX0X
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
/*	if( s->length()<8 ){
		ret.count=0;
		ret.score=0;
		return;
	}
	const size_t Lk = s->length()-8;
	const char * const ptr       = s->get_sequence();
	//put all kmers of template sequence into hash
	Simple_hash<int> h(Lk);
	for( size_t j=0; j<=Lk; ++j ){ 
		kmer[0] = ptr[j];
		kmer[1] = ptr[j+1];
		kmer[2] = ptr[j+3];
		kmer[3] = ptr[j+4];
		kmer[4] = ptr[j+5];
		kmer[5] = ptr[j+7];
		h.put(kmer, k, j);
	}

	//sum up max-kmer-score per position in query
	Kmer_seed_list::Iterator it = seed_list.begin();
	const Kmer_seed_list::Iterator end = seed_list.end();
	ret.count   = 0;
	ret.score  = 0.0f;
	float max_kmer_score = 0.0f;
	int last_pos = 0;
	if( it!=end ) last_pos = (*it)->position; 
	
	char mer[k];
	while( it!=end ){
		if( last_pos != (*it)->position ){
			ret.score += max_kmer_score;
			if(max_kmer_score>0.0f) ++ret.count;
			max_kmer_score = 0.0f;
			last_pos = (*it)->position;	
		}
		index2kmer((*it)->index, mer);
		if( h.exists(mer, k) ){
			max_kmer_score = std::max(max_kmer_score, (*it)->score);
		}
		++it;
	}
	ret.score += max_kmer_score;
	if(max_kmer_score>0.0f) ++ret.count;
}

//equivalent to match_full method of base_table
void Kmer::get_full_kmer_score_with(Sequence *s, Kmer::kmer_score &ret){
	if( s->length()<k ){
		ret.count=0;
		ret.score=0;
		return;
	}
	const size_t Lk = s->length()-k;
	const char * const ptr       = s->get_sequence();
	//put all kmers of template sequence into hash
	Simple_hash<int> h(Lk);
	for( size_t j=0; j<=Lk; ++j ){ 
		//for(size_t i=j; i<j+k; ++i) std::cout << INT2AA[i];
		//std::cout << " " << j << std::endl;
		h.put(&ptr[j], k, j);
	}

	//sum up max-kmer-score per position in query
	Kmer_seed_list::Iterator it = seed_list.begin();
	const Kmer_seed_list::Iterator end = seed_list.end();
	ret.count = 0;
	ret.score = 0.0f;
	
	char mer[k];
	while( it!=end ){
		index2kmer((*it)->index, mer);
		if( h.exists(mer, k) ){
			//std::cerr << "Match: " << (*it)->score << " ";
			//print_kmer( (*it)->index, std::cerr );
			//std::cerr << std::endl;
			ret.score += (*it)->score;
			++ret.count;
		}
		++it;
	}
}*/

size_t Kmer::count_idents(Sequence *q, Sequence *s, size_t k){
	size_t ret=0;
	if( (q->length()<k) || (s->length()<k) ) return 0;
	const size_t Lq = q->length()-k;
	const size_t Ls = s->length()-k;
	const char * const qseq      = q->get_sequence();
	const char * const sseq      = s->get_sequence();
	//build hash
	Simple_hash< std::list<int>* > h(Lq);
	for( size_t j=0; j<=Lq; ++j ){ 
		if( h.exists( &qseq[j], k ) ){
			std::list<int> *l = *(h.get(&qseq[j], k));
			l->push_back(j);
		}else{ 
			std::list<int> *l = new std::list<int>;
			l->push_back(j);
			h.put(&qseq[j], k, l); 
		}
	}
	for( size_t j=0; j<=Ls; ++j )
		if( h.exists( &sseq[j], k ) ){
			std::list<int> *l = *(h.get(&sseq[j], k));
			std::list<int>::iterator it  = l->begin();
			std::list<int>::iterator end = l->end();
			while( it!=end ){
				++it;
				++ret;
			}
		}
	Simple_hash< std::list<int>* >::Iterator it = h.begin();
	Simple_hash< std::list<int>* >::Iterator end = h.end();
	while( it!=end ){
		delete it.getValue();
		++it;
	}
	return ret;
}

std::ostream& Kmer::print_idents(Sequence *q, Sequence *s, size_t k, std::ostream &out){
	if( (q->length()<k) || (s->length()<k) ) return out;
	const size_t Lq = q->length()-k;
	const size_t Ls = s->length()-k;
	const char * const qseq      = q->get_sequence();
	const char * const sseq      = s->get_sequence();
	//build hash
	Simple_hash< std::list<int>* > h(Lq);
	for( size_t j=0; j<=Lq; ++j ){ 
		if( h.exists( &qseq[j], k ) ){
			std::list<int> *l = *(h.get(&qseq[j], k));
			l->push_back(j);
		}else{ 
			std::list<int> *l = new std::list<int>;
			l->push_back(j);
			h.put(&qseq[j], k, l); 
		}
	}
	for( size_t j=0; j<=Ls; ++j )
		if( h.exists( &sseq[j], k ) ){
			std::list<int> *l = *(h.get(&sseq[j], k));
			std::list<int>::iterator it  = l->begin();
			std::list<int>::iterator end = l->end();
			while( it!=end ){
				out << j << "  " << (*it) << std::endl;
				++it;
			}
		}
	Simple_hash< std::list<int>* >::Iterator it = h.begin();
	Simple_hash< std::list<int>* >::Iterator end = h.end();
	while( it!=end ){
		delete it.getValue();
		++it;
	}
	return out;
}



size_t Kmer::get_best_diag(size_t *diags, const size_t len){
	const int W = 100;
	size_t s_min=0, s_plus=0, s_tri=0, s_tri_max;
	//holds the diag with the maximum score
	size_t d_max;
	//init s_tri with the score at position W:
	// left(minus) window: 0..W-1 W W+1..W+(W-1)+1 right(plus) window
	for(int j=0; j<=W; ++j){
		s_tri += (j+1) * diags[j]; 
		if( j<W )
			s_tri += (W-j) * diags[W+j+1] ; 
	}
	for(int j=0; j<=W; ++j){
		s_min  += diags[j];
		s_plus += diags[W+j+1];
	}
	//printf("STRI:%i SMIN:%i SPLU:%i\n", s_tri, s_min, s_plus);
	s_tri_max = s_tri;
	d_max = W;
	for( size_t i=101; i<len-W; ++i){
		s_tri += -s_min +s_plus;
		if( s_tri>s_tri_max ){
			s_tri_max = s_tri;
			d_max     = i;
		}
		s_min  += diags[i] - diags[i-W-1];
		s_plus += -diags[i] + diags[i+W+1];
	}
	//printf("MAXIMUM DIAGONAL: %i with score %i\n", d_max, s_tri_max);
	return d_max;
}

void Kmer::_print_diags(size_t *diags, const size_t len){
	for( size_t i=0; i<len; ++i)printf("%Zu %Zu\n", i, diags[i]);
}


size_t Kmer::create_kmer_list_fast(Sequence *q) throw (std::exception){
	if (m_4mer != NULL) return create_kmer_list_from_matrix(q);
	
	int kmer[k];
	float smax[k];
	float corr[k];
	
	seed_list.reset();
	init_to_sequence(q);
	if( q->length()<k ) return 0;
	const size_t Lk        = q->length()-k;
	
	kmer_count             = 0;	
//	std::cout << q->get_header() << "\n";

//	std::cout << "|| " << q->index() << "||----------------------------------------------------------------------------------------------------\n";

	for( size_t i=0; i<=Lk; ++i ){
		q->get_kmer_at(kmer, k, i);
		memcpy(&corr, &score_corr[i], k*sizeof(float));
		float sum = 0.0;
		for (int z = 0; z < k; z++){
			sum += corr[z];
		}
//		int kmer_count_old = kmer_count;
		switch(k){
			case 4:{
				_init_smax(kmer, smax, corr);
				_add_similar_4mers_fast(i, kmer[0], kmer[1], kmer[2], kmer[3], smax, corr);
				break;
			}
			case 6:{
				_init_smax(kmer, smax, corr);
				_add_similar_6mers_fast(i, kmer[0], kmer[1], kmer[2], kmer[3], kmer[4], kmer[5], smax, corr);
				break;
			}
			default:{
				throw MyException("Fast k-mer generation not implemented for k=%zu!",k);
			}
		}
	}
	return kmer_count;	
}



//PATTERNS: k=6  XX0XXX0X
//				k=4  XX0XX
size_t Kmer::create_kmer_list_fast_spaced(Sequence *q) throw (std::exception){
	if (m_4mer != NULL) return create_kmer_list_spaced_from_matrix(q);

	int kmer[k];
	float smax[k];
	float corr[k];
	
	init_to_sequence(q);
	seed_list.reset();
	if( q->length()<8 ) return 0;
	const size_t Lk         = q->length()-8;
	const char * const ptr  = q->get_sequence();
	kmer_count       = 0;
	for( size_t i=0; i<=Lk; ++i ){
		switch(k){
			case 4:{
				corr[0] = score_corr[i];
				corr[1] = score_corr[i+1];
				corr[2] = score_corr[i+3];
				corr[3] = score_corr[i+4];

				kmer[0] = q->get_kmer_entry_at(i);
				kmer[1] = q->get_kmer_entry_at(i+1);
				kmer[2] = q->get_kmer_entry_at(i+3);
				kmer[3] = q->get_kmer_entry_at(i+4);
				_init_smax(kmer, smax, corr);
				_add_similar_4mers_fast(i, kmer[0], kmer[1], kmer[2], kmer[3], smax, corr);
				break;
			}
			case 6:{
				corr[0] = score_corr[i];
				corr[1] = score_corr[i+1];
				corr[2] = score_corr[i+3];
				corr[3] = score_corr[i+4];
				corr[4] = score_corr[i+5];
				corr[5] = score_corr[i+7];

				kmer[0] = q->get_kmer_entry_at(i);
				kmer[1] = q->get_kmer_entry_at(i+1);
				kmer[2] = q->get_kmer_entry_at(i+3);
				kmer[3] = q->get_kmer_entry_at(i+4);
				kmer[4] = q->get_kmer_entry_at(i+5);
				kmer[5] = q->get_kmer_entry_at(i+7);

				_init_smax(kmer, smax, corr);
				_add_similar_6mers_fast(i, kmer[0], kmer[1], kmer[2], kmer[3], kmer[4], kmer[5], smax, corr);
				break;
			}
			default:{
				throw MyException("Fast spaced k-mer generation implemented for k=6 or k=4 only!");
			}
		}
	}
	return kmer_count;	
}


void Kmer::_add_similar_6mers_fast(const size_t pos, int kmer0, int kmer1, int kmer2, int kmer3, int kmer4, int kmer5, float* smax, float* corr){
	register size_t index0, index1, index2;
	register float score0, score1, score2;
	register size_t i,j,l;
	const size_t pow2 = aa_powers[2];
	const size_t pow4 = aa_powers[4];
	const int d1 = kmer0 + kmer1*AMINOACID_DIM;
	const int d2 = kmer2 + kmer3*AMINOACID_DIM;
	const int d3 = kmer4 + kmer5*AMINOACID_DIM;
	for( i=0; i<pow2; ++i ){	
		score0 = _dimers[d1][i].score + corr[0] + corr[1];
		index0 = _dimers[d1][i].idx ;
		if( score0<smax[1] ) break;
		for( j=0; j<pow2; ++j ){
			score1 = score0 + _dimers[d2][j].score + corr[2] + corr[3];
			index1 = index0 + _dimers[d2][j].idx*pow2;
			if( score1<smax[3] ) break;
			for( l=0; l<pow2; ++l ){
				score2 = score1 + _dimers[d3][l].score  + corr[4] + corr[5];
				index2 = index1 + _dimers[d3][l].idx*pow4;
				if( score2>k_thr ){
					seed_list.append(index2, score2, pos);
					++kmer_count;
				}else break;
			}
		}
	}
}

void Kmer::_add_similar_4mers_fast(const size_t pos, int kmer0, int kmer1, int kmer2, int kmer3, float* smax, float* corr){
	register size_t index0, index1;
	register float score0, score1;
	register size_t i,j;
	const size_t pow2 = aa_powers[2];	
	const int d1 = kmer0 + kmer1*AMINOACID_DIM;
	const int d2 = kmer2 + kmer3*AMINOACID_DIM;
	
	for( i=0; i<pow2; ++i ){	
		score0 = _dimers[d1][i].score + corr[0] + corr[1];
		index0 = _dimers[d1][i].idx ;
		if( score0<smax[1] ) break;
		for( j=0; j<pow2; ++j ){
			score1 = score0 + _dimers[d2][j].score + corr[2] + corr[3];
			index1 = index0 + _dimers[d2][j].idx*pow2;
			if( score1>k_thr ){
				seed_list.append(index1, score1, pos);
				++kmer_count;
			}else break;
		}
	}
}
size_t Kmer::create_kmer_list_from_matrix(Sequence *q) throw (std::exception){
	int kmer[k];
	float corr[k];
	
	kmer_count = 0;
	seed_list.reset();
	if( q->length()<k ) return 0;
	const size_t Lk        = q->length()-k;
		
	init_to_sequence(q);
	for( size_t i=0; i<=Lk; ++i ){
		q->get_kmer_at(kmer, k, i);
		memcpy(&corr, &score_corr[i], k*sizeof(float));
		switch(k){
			case 4:{
				_add_similar_4mers_from_matrix(i, kmer, corr);
				break;
			}
			case 6:{
				_add_similar_6mers_from_matrix(i, kmer, corr);
				break;
			}
			default:{
				throw MyException("k-mer generation from k-mer similarity matrix not implemented for k=%zu!",k);
			}
		}
	}
//	std::cout << "kmer count for sequence: " << kmer_count << "\n_________________________________________\n\n";
	return kmer_count;
}

size_t Kmer::create_kmer_list_spaced_from_matrix(Sequence* q) throw (std::exception){
	int kmer[k];
	float smax[k];
	float corr[k];

	init_to_sequence(q);
	seed_list.reset();
	if( q->length()<8 ) return 0;
	const size_t Lk         = q->length()-8;
	const char * const ptr  = q->get_sequence();
	kmer_count       = 0;

	for( size_t i=0; i<=Lk; ++i ){
		switch(k){
			case 4:{
				corr[0] = score_corr[i];
				corr[1] = score_corr[i+1];
				corr[2] = score_corr[i+3];
				corr[3] = score_corr[i+4];

				kmer[0] = q->get_kmer_entry_at(i);
				kmer[1] = q->get_kmer_entry_at(i+1);
				kmer[2] = q->get_kmer_entry_at(i+3);
				kmer[3] = q->get_kmer_entry_at(i+4);

				_add_similar_4mers_from_matrix(i, kmer, corr);
				break;
			}
			case 6:{
				corr[0] = score_corr[i];
				corr[1] = score_corr[i+1];
				corr[2] = score_corr[i+3];
				corr[3] = score_corr[i+4];
				corr[4] = score_corr[i+5];
				corr[5] = score_corr[i+7];

				kmer[0] = q->get_kmer_entry_at(i);
				kmer[1] = q->get_kmer_entry_at(i+1);
				kmer[2] = q->get_kmer_entry_at(i+3);
				kmer[3] = q->get_kmer_entry_at(i+4);
				kmer[4] = q->get_kmer_entry_at(i+5);
				kmer[5] = q->get_kmer_entry_at(i+7);

				_add_similar_6mers_from_matrix(i, kmer, corr);
				break;
			}
			default:{
				throw MyException("Fast spaced k-mer generation implemented for k=6 or k=4 only!");
			}
		}
	}
	return kmer_count;
}

void Kmer::_add_similar_6mers_from_matrix(const size_t seq_pos, int* kmer, float* corr){
	
	int kmer_idx = kmer2index(kmer, 0, 6);
	int _4mer_idx = kmer2index(kmer, 0, 4);
	int _2mer_idx = kmer2index(kmer, 4, 6);
	
	Simple_list<Kmer_matrix::entry_score>* l4 = m_4mer->get_kmer_list(_4mer_idx);
	if (l4->length() == 0) return;
	
	Simple_list<Kmer_matrix::entry_score>::Iterator it4 = l4->begin();
	Simple_list<Kmer_matrix::entry_score>::Iterator end4 = l4->end();
	
	Simple_list<Kmer_matrix::entry_score>* l2 = m_2mer->get_kmer_list(_2mer_idx);
	if (l2->length() == 0) return;
	
	Simple_list<Kmer_matrix::entry_score>::Iterator begin2 = l2->begin();
	Simple_list<Kmer_matrix::entry_score>::Iterator end2 = l2->end();
	
	float _4corr = corr[0] + corr[1] + corr[2] + corr[3];
	float _2corr = corr[4] + corr[5];
	
	int _4thr = k_thr - ((*begin2).score + _2corr);
	int sim_6mer_index;
	float score;
	
	while(it4 != end4) {
		if (((*it4).score + _4corr) < _4thr) break;
		
		for (Simple_list<Kmer_matrix::entry_score>::Iterator it2 = begin2; it2 != end2; it2++){
			score = (*it4).score + (*it2).score + _4corr + _2corr;
			if ( score < k_thr) break;
			
			sim_6mer_index = (*it4).kmer_idx + (*it2).kmer_idx * aa_powers[4];
			seed_list.append(sim_6mer_index, score, seq_pos);
			++kmer_count;
		}
		++it4;
	}
}

void Kmer::_add_similar_4mers_from_matrix(const size_t seq_pos, int* kmer, float* corr){
	int kmer_idx = kmer2index(kmer, 0, 4);
//	std::cout << kmer_idx << "\n";
	Simple_list<Kmer_matrix::entry_score>* l = m_4mer->get_kmer_list(kmer_idx);
	if (l->length() == 0) return;
	
	Simple_list<Kmer_matrix::entry_score>::Iterator it = l->begin();
	Simple_list<Kmer_matrix::entry_score>::Iterator end = l->end();
	
	float _4corr = corr[0] + corr[1] + corr[2] + corr[3];
	
	while( it!=end) {
		if (((*it).score + _4corr) < k_thr) break;
		seed_list.append((*it).kmer_idx, ((*it).score + _4corr), seq_pos);
		++kmer_count;
		++it;
	}
}


int Kmer::create_idents_list(Sequence *q){
	int kmer[k];
	init_to_sequence(q);
	seed_list.reset();
	if( q->length()<k ) return 0;
	const size_t Lk        = q->length()-k;
	for( size_t i=0; i<=Lk; ++i ){
		q->get_kmer_at(kmer, k, i);
		seed_list.append(kmer2index(kmer, 0, k), 1.0, i);
	}
	return Lk;
}

/*void Kmer::index2kmer(size_t index, char *kmer){
	for(size_t i=0; i<k; ++i){
		size_t aa = index%21;
		kmer[i] = aa;
		index = index-aa;
		index /= AMINOACID_DIM;
	}
}*/

std::string Kmer::idx2kmer (int idx){
    std::string ret = "";
    int pos;
    for (int i = 0; i < k; i++){
        pos = idx % AMINOACID_DIM;
        ret += int2aa[pos];
        idx /= AMINOACID_DIM;
    }
    return ret;
}

size_t Kmer::kmer2index( const int *kmer, int begin, int end)const{
	size_t ret=0;
	for( size_t i=begin; i<end; ++i ) ret += kmer[i]*aa_powers[i-begin];
	return ret;
}

void Kmer::_print_kmer(char *ptr){
	for(size_t i=0; i<k; ++i)printf("%c", int2aa[(int)ptr[i]]);
}

std::ostream& Kmer::print_kmer(size_t index, std::ostream &out){
	for(size_t i=0; i<k; ++i){
		size_t aa = index%21;
		out << int2aa[aa];
		index = index-aa;
		index /= AMINOACID_DIM;
	}
	return out;
}
