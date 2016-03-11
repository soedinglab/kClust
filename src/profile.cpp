#include "profile.h"


/*Profile::Profile(	const char* const header, 
			const char* const sequence, 
			size_t index
			) throw (std::exception): 
			Sequence (header, 
			sequence, 
			index){
	members = new std::list<Member>();
}

Profile::Profile(char* header, 
				char* sequence, 
				const int* const aa2int, 
				const char* const int2aa, 
				size_t ADIM, 
				size_t index
				) throw (std::exception): 
				Sequence(header,
					sequence,
					aa2int,
					int2aa,
					ADIM,
					index){
	members = new std::list<Member>();
}*/
					
Profile::Profile(	const char* const header, 
			const char* const sequence, 
			size_t index,
			bool sc,
			Matrix *m, 
			bool dummy) throw (std::exception):
			Sequence(header, sequence, sc, index, m, dummy){
	members = new std::list<int>();
	profile = NULL;
	max_scores = NULL;
}

Profile::~Profile(){
//	std::cout << "Deleting profile ";
	if (profile != NULL) {
//		std::cout << " and matrix";
		for (size_t i = 0; i < len; i++){
				delete [] profile[i];
		}
		delete [] profile;
	delete [] max_scores;
	}
//	std::cout << "\n";
	delete members;
}

void Profile::init_profile_matrix(float** p){
	profile = p;
	max_scores = new float[len];
	for (int i = 0; i < len; i++){
		float max = 0.0f;
		for (int j = 0; j < m->get_aa_dim(); j++){
			if (profile[i][j] > max) max = profile[i][j];
		}
		max_scores[i] = max;
	}
}


void Profile::add_member(int orig_index){
	members->push_back(orig_index);
}

float Profile::get_score_with (int seq_pos, char aa) { 
	return profile[seq_pos][aa];
}

size_t Profile::get_memory_usage(){
	
//	std::cout << "Profile memory usage";
	
	size_t ret = Sequence::get_memory_usage();
	if (profile != NULL){
//		std::cout << " with matrix";
	    // profile matrix
	    ret += len * m->get_aa_dim() * sizeof(float);
	    // pointer to the matrix
	    ret += sizeof(size_t);
	    // array of second level pointers to the matrix
	    ret += len * sizeof(size_t);
	    // max_scores array
	    ret += len * sizeof(float);
	}
	// members array
/*	for (std::list<Member>::iterator it = members->begin(); it != members->end(); it++){
	    // Member object memory usage
	    ret += (*it).get_memory_usage();
	    // pointer to the next element
	    ret += sizeof(size_t);
	}*/
	ret += members->size() * sizeof(size_t);
//	std::cout << ": "<< ret << "\n";
	return ret;
}

void Profile::clear(){
    if (profile != NULL) {
		for (size_t i = 0; i < len; i++){
		    delete [] profile[i];
		}
		delete [] profile;
		delete [] max_scores;
    }
    profile = NULL;
    max_scores = NULL;
}

void Profile::calc_score_correction (int pseudocounts){
	correction = new float [len];
	memset(correction, 0.0f, len*sizeof(float));
}

float* Profile::get_score_correction(){
	return correction;
}
