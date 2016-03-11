#ifndef CM_SEQUENCE_PROFILE_H
#define CM_SEQUENCE_PROFILE_H

#include <list>
#include "sequence.h"

class Profile : public Sequence {
	public:
	
		Profile(	const char* const header, 
				const char* const sequence, 
				size_t index,
				bool sc,
				Matrix *m,
				bool dummy
			) throw (std::exception);
			
		virtual ~Profile();
		
/*		class Member{
			public:
				Member(int orig_idx, char* hdr): orig_index(orig_idx), header(hdr){}
				Member(const Member &m){
					orig_index = m.orig_index;
					header = new char [strlen(m.header)+1];
					strcpy(header, m.header);
				}
				~Member(){ 
					delete[] header;
				}
				bool operator==(const Member &m){ 
					return m.orig_index == orig_index;
				}
				Member& operator=(const Member &m){ 
					orig_index = m.orig_index;
					header = new char [strlen(m.header)+1];
					strcpy(header, m.header);
					return *this;
				}
				bool operator<(const Member &m){ 
					return orig_index < m.orig_index;
				}
				bool operator>(const Member &m){ 
					return orig_index > m.orig_index;
				}
				size_t get_memory_usage(){
				    // orig_index
				    size_t ret = sizeof(int);
				    // header
				    ret += strlen(header) * sizeof(char);
				}
			
			public:
				int orig_index;
				char* header;
		};*/
		
		void init_profile_matrix(float** profile);
		
		void add_member(int orig_index);

		std::list<int>* get_members(){ return members; };
		
		virtual inline const float** const get_scoring_matrix() { return (const float**) profile;}
		virtual inline const float * const get_max_scores() { return max_scores; }
		virtual inline void get_kmer_at (int * kmer, int k, int pos) { 
			for (int i = 0; i < k; i++) *(kmer++) = pos+i; 
		}
		virtual inline int get_kmer_entry_at (int pos) { return pos; }
		virtual inline int get_aa_entry_at (int pos) { return _seq[pos]; }
		virtual float get_score_with (int seq_pos, char aa);
		// clears profile data when saved as representative
		virtual void clear();
		virtual size_t get_memory_usage();
		virtual float* get_score_correction();
	
	private:
		virtual void calc_score_correction (int pseudocounts = 20);
		// profile has the size len x AMINOACID_DIM
		float** profile;
		float* max_scores;
		// members of the cluster from which the profile was constructed
		std::list<int> * members;
};

#endif
