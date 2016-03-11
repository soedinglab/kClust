#include "fasta_db_reader.h"
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
/* TODO die Buffer Groessen genauer anpassen (während dem ersten durchlauf den speicherbedarf messen?)
 */
Fasta_db_reader::Fasta_db_reader( std::string f, Matrix *matrix, bool _sc, bool profile_query, std::string adir, std::string c_old_file, std::string h_old_file):
	file(f),
	_profile_query(profile_query),
	m(matrix),
	sc(_sc),
	AMINOACID_DIM( matrix->get_aa_dim()),
	alndir(adir),
	clusters_old_file(c_old_file),
	headers_old_file(h_old_file)
	{
				
	hbuffer = new char[READ_HBUFFER_SIZE];
	sbuffer = new char[READ_SBUFFER_SIZE];
	if (_profile_query){
		h_old_buffer = new char[READ_HBUFFER_SIZE];
		pbuffer = new char[READ_PBUFFER_SIZE];
	}

	max_seq_len   = 0;
	min_seq_len   = INT_MAX;
	seq_count     = 0;	
	orig_db_size  = 0;
	_seeks_       = 0;
	_hseeks_      = 0;
	clusters_old  = 0;	
	current_idx   = 1;

	chars_d              = 0.0;
	sum_log_template_len = 0.0;

	srand( time(NULL) );

	in.open(file.c_str());
	if( in.fail() ) throw MyException("Cannot read '%s'", file.c_str());
	
	if (profile_query){
		if (headers_old_file == "") std::cerr<<"No path to the old headers file for profile query passed to the database reader!";
		in_h.open(headers_old_file.c_str());
		if( in_h.fail() ) throw MyException("Cannot read '%s'", file.c_str());
	}

	// move to first header line 
	while( in.good() ){
		char c = in.get();
		if( c=='>' ){ in.unget(); break; }	
		if( isspace(c) ) continue;
		throw MyException("File '%s' is not in correct FastA format!", file.c_str());
	}
	if (profile_query && alndir == "") std::cerr<<"No input directory for profile query passed to the database reader!";
	if (profile_query && clusters_old_file == "") std::cerr<<"No path to the old clusters file for profile query passed to the database reader!";
	
	_start_of_data = in.tellg();
	_start_of_data_h = in_h.tellg();
	_init_limits_();
	_init_seeks_();
	reset();
}

Fasta_db_reader::~Fasta_db_reader(){
	delete [] _seeks_;
	delete [] hbuffer;
	delete [] sbuffer;
	delete [] idx2orig;
	delete [] orig2idx;
	if (_profile_query){
		delete [] h_old_buffer;
		delete [] clusters_old;
		delete [] _hseeks_;
		delete [] pbuffer;
	}
}

void Fasta_db_reader::reset() throw (std::exception){
	in.clear();
	in.seekg( _start_of_data, std::ios::beg );
	current_idx = 1;
}

void Fasta_db_reader::reset_headers() throw (std::exception){
	in_h.clear();
	in_h.seekg( _start_of_data_h, std::ios::beg );
}

void Fasta_db_reader::_init_limits_() throw (std::exception){

	while(!in.eof()){
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		size_t hb = in.gcount();
		
		if( hb==(READ_HBUFFER_SIZE-1) ) 
			throw MyException("Error reading '%s'\nLong header line found, increase READ_HBUFFER_SIZE!", file.c_str());	
			
		if( hbuffer[0] == '>' ){
			in.get(sbuffer, READ_SBUFFER_SIZE, '>');
			size_t sb = in.gcount();
			if( sb==(READ_SBUFFER_SIZE-1) )
				throw MyException("Error reading '%s'\nLong sequence found, increase READ_SBUFFER_SIZE!", file.c_str());
			++seq_count;
			const size_t slen = Sequence::get_length( sbuffer, sb);
			if( slen > max_seq_len ) max_seq_len = slen;
			if( slen < min_seq_len ) min_seq_len = slen;
			chars_d +=(double)slen;
			sum_log_template_len += log(slen);
		}else throw MyException("File '%s' is not in correct FastA format!", file.c_str());
	}
	if(seq_count==0) throw MyException("File '%s' contains no sequences in FastA format!", file.c_str());
	std::cout << "DB size: " << seq_count << "\n";
	// for profiles, seq_count is the number of the consensus sequences, orig_db_size the overall number of sequences to cluster
	// for non profile queries, seq_count=orig_db_size

	if (! _profile_query) {
		orig_db_size = seq_count;

		idx2orig = new size_t [seq_count+1];
		for (int i = 0; i < (seq_count+1); i++) idx2orig[i] = i;
		
		orig2idx = new size_t [seq_count+1];
		for (int i = 0; i < (seq_count+1); i++) orig2idx[i] = i;
	}
	else _init_profile_settings_();
	
//	std::cout << "seq_count: " << seq_count << "\n";
//	std::cout << "orig_db_size: " << orig_db_size << "\n";
		
	//std::cerr << "Number of sequences in database: " << seq_count   << std::endl;
	//std::cerr << "Minimum sequence length        : " << min_seq_len << std::endl;
	//std::cerr << "Maximum sequence length        : " << max_seq_len << std::endl;
}

// Counts the real size of the database for the profile query and fills orig2idx array
void Fasta_db_reader::_init_profile_settings_(){
	
	std::ifstream in_c(clusters_old_file.c_str());
	
	in_c.getline(hbuffer, READ_HBUFFER_SIZE);
	// read the database size
	hbuffer[0] = ' ';
	orig_db_size = atoi(hbuffer);
	// fill clusters_old array
	clusters_old = new std::list<size_t> [orig_db_size+1];
	
//	std::cout << "Init old clusters\n";

	while (!in_c.eof()){
		in_c.getline(hbuffer, READ_HBUFFER_SIZE);
		if (in_c.gcount() == 0) continue;
		// orig_idx of sequence and representant of its cluster
		size_t seq = atoi(strtok(hbuffer," "));
		size_t rep = atoi(strtok(NULL," "));
		clusters_old[rep].push_back(seq);
	}
	in_c.close();
	
//	std::cout << "\nInit idx2orig, orig2idx\n";

	idx2orig = new size_t [seq_count+1];
	orig2idx = new size_t [orig_db_size+1];
	
	// fill idx2orig array
	in.clear();
	in.seekg( _start_of_data, std::ios::beg );
	
	int idx = 1;
	int orig_idx = 1;
	int line_pos = 1;
	
	idx2orig[0] = 0;
	for (int i = 0; i < (orig_db_size+1); i++) {orig2idx[i] = 0;}
	
	while(!in.eof()){
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		
		hbuffer[0] = ' ';
		orig_idx = atoi(hbuffer);
		
//		std::cout << hbuffer << "\n";

		if (idx>seq_count) std::cout << "warning: idx = "<<idx <<"\n";
		if (orig_idx>orig_db_size) std::cout << "warning: orig_idx = "<<orig_idx <<"\n";
		
		idx2orig[idx] = orig_idx;
		orig2idx[orig_idx] = idx;
		
//		std::cout << orig_idx << " " << idx << "\n";

		idx++;
	}
} 


/*std::ostream& Fasta_db_reader::analyze(std::ostream &out) throw (std::exception){
	reset();
	unsigned long counts[AMINOACID_DIM];
	for( size_t i=0; i<AMINOACID_DIM; ++i){
		counts[i] = 0L;
	}
	size_t seqs       =  0;
	size_t min_length =  INT_MAX;
	size_t max_length =  0;
	int *lengths      = new int[seq_count];
	long sum          = 0L;
	double sum_sqrt   = 0.0f;
	Sequence *s=0;
	while( !in.eof() ){
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		s = new Sequence(hbuffer, sbuffer, 0, m, false);
		if( (s->length()) > max_length ) max_length = s->length();
		if( (s->length()) < min_length ) min_length = s->length();
		lengths[seqs] = s->length();
		++seqs;
		sum += s->length();
		sum_sqrt += s->length()*s->length();
		for( size_t i=0; i<s->length(); ++i){
			++counts[(int)s->get_sequence()[i]];
		}
		delete s;
	}

	std::sort(lengths, lengths+seqs);

	out << "=============================================" << std::endl;
	out << "Minimum sequence length   : " << min_length       << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Maximum sequence length   : " << max_length       << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Number of sequences       : " << seqs             << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Characters                : " << sum              << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Average sequence length   : " << sum/(double)seqs << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Average sqrt of length^2  : " << sqrt(sum_sqrt/(double)seqs) << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Median sequence length    : " << lengths[seqs/2]  << std::endl;
	out << "---------------------------------------------" << std::endl;
	out << "Amino acid counts, frequencies:"               << std::endl<< std::endl;
	
	char buffer[512];
	const float * const p_back = m->get_p_back();
	const char*  const int2aa = m->get_int2aa();
	for(uint i=0; i<AMINOACID_DIM; ++i){

		sprintf(buffer, "%c %-10li %-5.2f %+5.2f\n",  int2aa[i], counts[i], 100.0f*counts[i]/(float)sum,  100.0f*((counts[i]/(float)sum) - p_back[i]) );
		out << buffer;
	}
	out << "=============================================" << std::endl;
	
	delete [] lengths;
	reset();
	return out;
}*/

std::ostream& Fasta_db_reader::analyze_bins(std::ostream &out, int width){
	reset();
	const size_t bin_count = (max_seq_len/width)+1;
	size_t *bins           = new size_t[bin_count];
	for( size_t i=0; i<bin_count; ++i) bins[i]=0;
	while( !in.eof() ){
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		Sequence *s = new Sequence(hbuffer, sbuffer, false, 0, m, false);
		++bins[s->length()/width];
		delete s;
	}
	out << "#Sequence length in bins (width=" << width << ", #bins=" << bin_count << ")" << std::endl;
	for( size_t i=0; i<bin_count; ++i){
		if( bins[i]!=0 )
			out << i << " " << bins[i] << std::endl;
	}
	delete [] bins;
	reset();
	return out;
}

Sequence* Fasta_db_reader::get_next( bool dummy, bool without_header, bool without_matrix){
	if (! _profile_query) return get_next_sequence (dummy, without_header);
	else return get_next_profile (dummy, without_header, without_matrix);
}

Sequence* Fasta_db_reader::get_next_sequence (bool dummy, bool without_header){
	Sequence *ret = 0;
	if(!in.eof()){
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		
		if (without_header) sprintf( hbuffer, ">%zu", current_idx );
		
//		std::cout << sbuffer;
		
		
//		if (!dummy) ret = new Sequence(hbuffer, sbuffer, m, current_idx);
//		else ret = new Sequence(hbuffer, sbuffer, current_idx);
		ret = new Sequence(hbuffer, sbuffer, sc, current_idx, m, dummy);
		
		current_idx++;
	}
	return ret;
}

Profile* Fasta_db_reader::get_next_profile(bool dummy, bool without_header, bool without_matrix){
	Profile *ret = 0;
	if(!in.eof()){	
		//read consensus sequence
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		
		size_t orig_idx = idx2orig[current_idx];

		if( orig_idx<1 || orig_idx > orig_db_size ) {
			std::cout << "WARNING get_next_profile orig_idx out of range " << orig_idx << "(db size " << orig_db_size << ")\n";
			std::cout << "Current idx " << current_idx << "\n";
			std::cout << "HEADER\n" << hbuffer << "\n";
		}

		if (without_header) sprintf( hbuffer, ">%zu", orig_idx );
		
		Profile *ret = 0;

		ret = new Profile(hbuffer, sbuffer, current_idx, sc, m, dummy);
		
		if (!without_matrix){
		    float** profile_m = get_profile_matrix(orig_idx, ret->length());
		    ret->init_profile_matrix(profile_m);
		}
		
  		for (std::list<size_t>::iterator it = clusters_old[orig_idx].begin() ; it != clusters_old[orig_idx].end(); it++ ){
  			ret->add_member(*it);
  	  		if (*it > orig_db_size || *it < 1){
  	  			std::cout << "WARNING get_next_profile orig_idx " <<  *it << "\n";
  	  			std::cout << "for representative orig_idx " << orig_idx << "\n";
  	  			std::cout << "rep HEADER\n" << hbuffer << "\n";
  	  		}
  		}
			
		current_idx++;
		
		return ret;
	}
	return ret;
}

Sequence* Fasta_db_reader::get( const size_t idx, bool dummy, bool without_header, bool without_matrix ) throw (std::exception) {
	if (! _profile_query) return get_sequence( idx, dummy, without_header );
	else return get_profile( idx, dummy, without_header, without_matrix );
}

Sequence* Fasta_db_reader::get_sequence( const size_t idx, bool dummy, bool without_header ) throw (std::exception) {
	if( idx<1 || idx > seq_count ) 
		throw MyException( "get_sequence: index of query '%zu' is out of range (<%zu)!", idx, seq_count );
	in.clear();
	in.seekg(_seeks_[idx-1], std::ios::beg);
	in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');		
	in.get(sbuffer, READ_SBUFFER_SIZE, '>');
	
	if (without_header) sprintf( hbuffer, ">%zu", idx );
	
//	current_idx = idx+1;
	
	return new Sequence(hbuffer, sbuffer, sc, idx, m, dummy);
}

Profile* Fasta_db_reader::get_profile( const size_t idx, bool dummy, bool without_header, bool without_matrix ) throw (std::exception) {
	if( idx<1 || idx > seq_count ) 
		throw MyException( "get_profile: index of query '%zu' is out of range (<%zu)!", idx, seq_count );
	// get the representative sequence
	in.clear();
	in.seekg(_seeks_[idx-1], std::ios::beg);

	in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');		
	in.get(sbuffer, READ_SBUFFER_SIZE, '>');
	
	if (without_header) sprintf( hbuffer, ">%zu", idx );
	
	Profile *ret = 0;
	size_t orig_idx = idx2orig[idx];
	
	if( orig_idx<1 || orig_idx > orig_db_size ) {
		std::cout << "WARNING get_profile orig_idx out of range " << orig_idx << "(db size " << orig_db_size << ")\n";
		std::cout << "Current idx " << idx << "\n";
		std::cout << "HEADER\n" << hbuffer << "\n";
	}

	ret = new Profile(hbuffer, sbuffer, sc, idx, m, dummy);		

	if (!without_matrix){
		float** profile_m = get_profile_matrix(orig_idx, ret->length());
		ret->init_profile_matrix(profile_m);
	}
	
  	for (std::list<size_t>::iterator it = clusters_old[orig_idx].begin() ; it != clusters_old[orig_idx].end(); it++ ){
  		ret->add_member(*it);
  		if (*it > orig_db_size || *it < 1){
  			std::cout << "WARNING get_profile orig_idx " <<  *it << "\n";
  			std::cout << "for representative orig_idx " << orig_idx << "\n";
  			std::cout << "rep HEADER\n" << hbuffer << "\n";
  		}
  	}
	
	return ret;

}

float** Fasta_db_reader::get_profile_matrix(int orig_idx, int seqlength){
//	std::cout << "get_profile_matrix for orig_idx = " << orig_idx << "\n";
	std::string hhmakelogf = "/tmp/hhmake.log";
	// calculate profile for the msa with hhmake
	char fname [10];
	fname[9] = '\0';
	Converter::get_file_name(fname, orig_idx);
	char dir [4] = {fname[0], fname[1], '/', '\0'};	
	
	std::string d_str (dir);
	std::string f_str (fname);	

/*	std::string hhmake_cmd = "/cluster/bioprogs/hh/hhmake -pca 1.0 -pcm 2 -Blosum65 -i " + alndir + d_str + f_str + ".a3m -o " + alndir + d_str + f_str + ".hhm -pca 1.4 -pcb 2.5 &> " + hhmakelogf;
	
	int sysret = system(hhmake_cmd.c_str());
	if (sysret != 0) throw MyException ("Cannot execute '%s'", hhmake_cmd.c_str());
	*/
	// read profile	
	std::ifstream in_p((alndir + d_str + f_str + ".hhm").c_str() );
	if (! in_p.is_open()) throw MyException( "Cannot open profile file " + alndir + d_str + f_str + ".hhm !");
	
	float ** ret = new float * [seqlength];
	for (int i = 0; i < seqlength; i++){
		ret[i] = new float[AMINOACID_DIM];
	}
	
	do {
		in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	}while( pbuffer[0] != '#');
	
	// skip next line (NULL)
	in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	
	// TODO jedes mal auslesen oder voraussetzen dass die aminosäuren alphabetisch geordnet sind und ein instanzvariablen-array anlegen???
	// read the order of aminoacids in the file
	char local_int2aa[20];
	int pos = 7;
	int i = 0;
	while (pos < 47){
		local_int2aa[i] = pbuffer[pos];
		pos = pos + 2;
		++i;
	}
	
	in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
	
	// init the profile with 0.0
	for (int i = 0; i < seqlength; i++){
		for (int j = 0; j < AMINOACID_DIM; j++){
			ret[i][j] = -1.0f;
		}
	}
	
	const int*  const aa2int = m->get_aa2int();
	const float * const p_back = m->get_p_back();
	
	
	int seq_pos = 0;
	while (pbuffer[0]!= '/' &&  pbuffer[1]!= '/'){
		int aa_num = 0;
		// insert a delimiter '\t' before the first profile entry
		pbuffer[6] = '\t';
		char * tok = strtok (pbuffer,"\t");
		while ((tok != NULL) && (aa_num < 20)){
			tok = strtok (NULL,"\t");
			// * entry: 0.0 probability
			if (tok[0] == '*') std::cerr<<"!!! 0 PROBABILITY FOR " << f_str << ".hhm AT " << seq_pos << "," << aa2int[local_int2aa[aa_num]]<<"\n";//ret[seq_pos][aa2int[local_int2aa[aa_num]]] = 0.0f;
			// 0 entry: 1.0 probability
			else if (tok[0] == '0') ret[seq_pos][aa2int[local_int2aa[aa_num]]] = m->_log2(1.0f / p_back[aa2int[local_int2aa[aa_num]]]) / m->get_scale();
			// integer number entry: 0.0 < probability < 1.0
			else{
				int entry = atoi(tok);
				float p = pow(2.0f, -(entry/1000.0f));
				ret[seq_pos][aa2int[local_int2aa[aa_num]]] = m->_log2( p / p_back[aa2int[local_int2aa[aa_num]]]) / m->get_scale();
			}
			aa_num++;
		}
		
		in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
		in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
		in_p.getline(pbuffer, READ_PBUFFER_SIZE, '\n');
			
		seq_pos++;
	}
	
	in_p.clear();
	in_p.close();
	
	// remove hhm file to save memory
//	system(("rm " + alndir + d_str + f_str + ".hhm").c_str());
//	system(("rm " + hhmakelogf).c_str());
	return ret;
}


// TODO das wird zur Zeit aber nicht verwendet
// wenn doch, muesste man das auf Profile erweitern
Packed_sequence* Fasta_db_reader::get_next_packed_without_header(){
	Packed_sequence *ret = 0;
	if(!in.eof()){	
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		sprintf( hbuffer, ">%zu", current_idx );
		ret = new Packed_sequence(hbuffer, sbuffer, current_idx++);
	}
	return ret;
}

size_t Fasta_db_reader::get_index_of_next_sequence(){ 
	return current_idx;
}

bool Fasta_db_reader::has_next() {
	return !in.eof();
}

bool Fasta_db_reader::move_forward_by(const size_t n){
	bool ret = false;
	size_t count=0;
	while( !in.eof() && (count<n) ){	
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		++count;
		current_idx++;
	}
	if( !in.eof() && (count==n) ) ret = true; 
	return ret;		
}

void Fasta_db_reader::move_to( const size_t idx ) throw (std::exception) {
	if( idx<1 || idx > seq_count ) 
		throw MyException( "move_to: index of query '%i' is out of range (<%i)!", idx, seq_count );
	std::cout << "Moving to index " << idx << ", orig_idx " << idx2orig[idx] << "\n";
	in.clear();
	in.seekg(_seeks_[idx-1], std::ios::beg);
	current_idx = idx;
}

void Fasta_db_reader::seek_sort( const std::string outfile ) throw (std::exception) {
	reset();
	SeekSortStruct *unsorted = new SeekSortStruct[seq_count];
	int i=0;
	int skipped = 0;
	while(!in.eof()){
		unsorted[i].pos = in.tellg();
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');		
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		unsorted[i].len   = Sequence::get_length( sbuffer, in.gcount() );
		unsorted[i].index = i+1;
		++i;	
	}
	SeekSortStruct **sorted = _merge_seek_sort( unsorted, seq_count );
	std::ofstream out(outfile.c_str());
	if( !out ) throw MyException("Cannot open '%s' for writing!", outfile.c_str());
	for( size_t i=0; i<seq_count; ++i){
		in.clear();
		in.seekg(sorted[i]->pos);
		in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');	
		size_t hbytes = in.gcount();	
		in.get(sbuffer, READ_SBUFFER_SIZE, '>');
		size_t sbytes = in.gcount();
		
		// check if it is a X* sequence
		bool accept = false;
		size_t spos = 0;
		while((sbuffer[spos] != '\0') && !accept){
			accept = (sbuffer[spos] != 'X' && sbuffer[spos] != 'x' && sbuffer[spos] != '\n' && sbuffer[spos] != ' ');
			++spos;
		}
		if (!accept){
			 // skip this sequence
//			 std::cerr << "Skipping sequence containing only X's: " << hbuffer << "\n";
			++skipped;
			continue;
		}
		
		out.write(hbuffer, hbytes-1);
		out << std::endl;
		out.write(sbuffer, sbytes);
		//the last sequence in the file may not have a terminal newline
		if( sbuffer[sbytes-1]!='\n' )
		out << std::endl;
	}
	out.close();
	delete [] sorted;
	delete [] unsorted;
//	std::cerr << "Skipped " << skipped << " sequences containing only X's.\n"; 
}


Fasta_db_reader::SeekSortStruct** Fasta_db_reader::_merge_seek_sort( SeekSortStruct *unsorted , size_t len ){
	if( len<2 ){
		SeekSortStruct **ret = new SeekSortStruct*[1];
		ret[0] = unsorted; 
		return ret;
	}
	size_t l=0;
	size_t h=0;
	size_t m=(len-1)/2;	
	SeekSortStruct **lower = _merge_seek_sort(unsorted, m+1 );
	SeekSortStruct **upper = _merge_seek_sort(unsorted+m+1, len-m-1);
	SeekSortStruct **tmp   = new SeekSortStruct*[len];
 	while (l+h<len){ 
    	if( l<=m && (h>=len-m-1 || lower[l]->len > upper[h]->len ) ){
			tmp[l+h] = lower[l];	
			++l;
		}else{ 
			tmp[l+h] = upper[h];
			++h;
		}
	}
	delete [] lower; 
	delete [] upper; 
  	return tmp;
}


std::ostream& Fasta_db_reader::write_random_partition( std::ostream &out, size_t sequences ){
        for(size_t i=0; i<sequences; ++i){
                size_t idx = int((rand()/(double)RAND_MAX) * (seq_count-1));
                in.clear();
                in.seekg(_seeks_[idx]);
                in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
                size_t hbytes = in.gcount();
                in.get(sbuffer, READ_SBUFFER_SIZE, '>');
                size_t sbytes = in.gcount();
                out.write(hbuffer, hbytes-1);
                out << std::endl;
                out.write(sbuffer, sbytes);
        }
        return out;
}


void Fasta_db_reader::_init_seeks_(){
	reset();
	delete [] _seeks_;
	_seeks_ = new std::streampos[seq_count];
	size_t i=0;
	while(!in.eof()){
		_seeks_[i++] = in.tellg();
        in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
        in.get(sbuffer, READ_SBUFFER_SIZE, '>');
	}
	if (_profile_query) _init_hseeks_();
}

void Fasta_db_reader::_init_hseeks_(){
	reset_headers();
	delete [] _hseeks_;
	_hseeks_ = new std::streampos[orig_db_size+1];

	while(!in_h.eof()){
		std::streampos pos = in_h.tellg();
        in_h.getline(h_old_buffer, READ_HBUFFER_SIZE, '\n');
        if (in_h.gcount() == 0) continue;
        size_t orig_idx = atoi(strtok(h_old_buffer," "));
        _hseeks_[orig_idx] = pos;
	}
}

char* Fasta_db_reader::get_header_for (int orig_idx){
	if( orig_idx<1 || orig_idx > orig_db_size ) 
		throw MyException( "get_header_for: index of query '%zu' is out of range (<%zu)!", orig_idx, orig_db_size );
	in_h.clear();
	in_h.seekg(_hseeks_[orig_idx], std::ios::beg);
	
	in_h.getline(h_old_buffer, READ_HBUFFER_SIZE, '\n');
	strtok(h_old_buffer," ");
	return strtok(NULL,"\0");
		
}

Sequence* Fasta_db_reader::get_rand(){
	size_t idx = int((rand()/(double)RAND_MAX) * (seq_count-1));
	in.clear();
	in.seekg(_seeks_[idx]);
	in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');	
	in.get(sbuffer, READ_SBUFFER_SIZE, '>');
//	return new Sequence(hbuffer, sbuffer, aa2int, int2aa, idx);
	return new Sequence(hbuffer, sbuffer, sc, idx, m, false);
}


const size_t Fasta_db_reader::get_memory_usage()const{
	size_t ret = sizeof(Fasta_db_reader);
	//seeks
	ret += seq_count*sizeof(std::streampos) + sizeof(size_t);
	//hseeks
	ret += (orig_db_size+1)*sizeof(std::streampos) + sizeof(size_t);
	//header buffer
	ret += sizeof(char)*READ_HBUFFER_SIZE   + sizeof(size_t);
	//sequence buffer
	ret += sizeof(char)*READ_SBUFFER_SIZE   + sizeof(size_t);
	//idx2orig, orig2idx
	ret += sizeof(size_t) *(seq_count+1) *2;
	if (_profile_query){
		//hseeks
		ret += (orig_db_size+1)*sizeof(std::streampos) + sizeof(size_t);
		//header old buffer
		ret += sizeof(char)*READ_HBUFFER_SIZE   + sizeof(size_t);
		//profile buffer
		ret += sizeof(char)*READ_PBUFFER_SIZE   + sizeof(size_t);
		//clusters_old list
		for (int i = 0; i < orig_db_size+1; i++){
			// not exact, without linking pointers
			ret += sizeof(char) *clusters_old[i].size();
		}
	}
	return ret;
}


