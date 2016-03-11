#include "kmer_matrix.h"

Kmer_matrix::Kmer_matrix(int k, float scale):
	k(k),
	scale(scale){
		AMINOACID_DIM = 21;
		size = pow(AMINOACID_DIM, k);
		scores = new Simple_list<entry_score> *[size];
		for (int i = 0; i < size; i++){
			scores[i] = new Simple_list<entry_score>();
		}
		
		aa_powers = new size_t[10];
		size_t pow = 1;
		for( size_t i=0; i<10; ++i ){
			aa_powers[i] = pow;
			pow *= AMINOACID_DIM;
		}
		
		int2aa = new char[AMINOACID_DIM];
	    int2aa[0] = 'A';
	    int2aa[1] = 'R';
	    int2aa[2] = 'N';
	    int2aa[3] = 'D';
	    int2aa[4] = 'C';
	    int2aa[5] = 'Q';
	    int2aa[6] = 'E';
	    int2aa[7] = 'G';
	    int2aa[8] = 'H';
	    int2aa[9] = 'I';
	    int2aa[10] = 'L';
	    int2aa[11] = 'K';
	    int2aa[12] = 'M';
	    int2aa[13] = 'F';
	    int2aa[14] = 'P';
	    int2aa[15] = 'S';
	    int2aa[16] = 'T';
	    int2aa[17] = 'W';
	    int2aa[18] = 'Y';
	    int2aa[19] = 'V';
	    int2aa[20] = 'X';
	
	    aa2int = new int['Z'+1];
	    for(size_t i=0; i<='Z'; ++i) aa2int[i]=-1;
	    for(size_t i=0; i<AMINOACID_DIM; ++i){
	            aa2int[int2aa[i]] = i;
	    }
}

Kmer_matrix::~Kmer_matrix(){
	for (int i = 0; i < size; i++){
		delete scores[i];
	}
	delete[] scores;
}

Simple_list<Kmer_matrix::entry_score>* Kmer_matrix::get_kmer_list(int kmer_idx){
	if (kmer_idx >= size) std::cerr << "idx: " << kmer_idx << ", size: " << size << "\n"; 
	return scores[kmer_idx];
}

void Kmer_matrix::fill_kmer_matrix(std::string file){
	struct stat stFileInfo;
	int s; 
	s = stat(file.c_str(),&stFileInfo);
  	if(s != 0) throw MyException("File " + file + " does not exist.\n");
  		
	std::ifstream in(file.c_str());

    char buf[1000];

    int idx1;
    int idx2;
    float score;

    int c = 0;
    while(!in.eof()){
        if (c%1000000 == 0) std::cerr << ".";
        in.getline(buf, 1000);
        if (buf[0] == '\0') continue;

        idx1 = atoi(strtok(buf,"\t"));
        idx2 = atoi(strtok(NULL,"\t"));
        score = atof(strtok(NULL,"\t")) / scale;
        
        Kmer_matrix::entry_score e (idx2, score);
        scores[idx1]->append(e);
        c++;
    }

    in.clear();
    in.close();
}