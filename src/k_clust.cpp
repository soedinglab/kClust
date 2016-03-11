/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include <boost/regex.hpp>
#include <dirent.h>
#include <sys/types.h> 
#include <fstream>
#include <iostream>

#include <unistd.h>
#include <ios>
#include <string>


#include "k_clust.h"


K_Clust::K_Clust( Params *p ) throw (std::exception):params(p), matrix( p->matrix_type ){
	node_count                     = 0;
	db_reader                      = 0;
	block_man                      = 0;
	clusters                       = 0;
	table                          = 0;
	_6mer                          = 0;
	_4mer                          = 0;
	_3mer                          = 0;
	_4meraln                       = 0;
	_4mer_table                    = 0;
	_4mer_matrix                   = 0;
	_2mer_matrix                   = 0;
	offset_memory_usage            = 0;

	kmer6_av = 0.0;
	kmer4_av = 0.0;
	seqs_6 = 0;
	seqs_4 = 0;
	matches_6mer = 0;
	matches_4mer = 0;

	benchmark_queue_size = 10000;

	q_pref = new std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare> ();
	q_kdp = new std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare> ();
	q_seq = new std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare> ();
}

K_Clust::~K_Clust(){
	delete db_reader;
	delete clusters;
	delete block_man;
	delete table;
	delete _6mer;
	delete _4mer;
	delete _3mer;
	delete _4mer_matrix;
	delete _2mer_matrix;
	delete _4mer_table;
	delete _4meraln;
	delete q_pref;
	delete q_kdp;
	delete q_seq;
}

void K_Clust::_init() throw (std::exception){

	std::cerr << "Sorting database...";
	Fasta_db_reader dbr( params->get_dbfile(), &matrix, false, false );
	dbr.seek_sort( params->get_dbsortedfile() );	
	std::cerr << "ok" << std::endl;
	std::cerr << "======================================================" << std::endl;

	// for the profile query the database file is only temporary and is removed
	/*    if (params->profile_query){
	      std::string rmcmd = "rm " + params->get_dbfile();
	      int ret = system(rmcmd.c_str());
	      if (ret != 0) throw MyException ("Cannot execute '%s'", rmcmd.c_str());
	      int ret = system(("cp " + params->get_dbsortedfile_old() + " " + params->get_working_dir()).c_str());
	      if (ret != 0) throw MyException ("Cannot execute '%s'", ("cp " + params->get_dbsortedfile_old() + " " + params->get_working_dir()).c_str());
	      }*/

	std::cerr << "Init..." << std::endl;

	std::cerr << "   Database reader...";
	db_reader   = new Fasta_db_reader( params->get_dbsortedfile(), &matrix, params->score_correction, params->profile_query, params->get_aln_dir_old(), params->get_clustersfile_old(), params->get_headersfile_old() );
	std::cerr << "ok" << std::endl;

	std::cerr << "   Clusters...";
	clusters    = new Clusters( db_reader->get_sequence_count(), &matrix );
	std::cerr << "ok" << std::endl;

	std::cerr << "   Block manager...";
	block_man   = new BlockManager( db_reader, clusters ); 
	std::cerr << "ok" << std::endl;

	std::cerr << "   Index-table for filtering...";
	table       = new Base_table( db_reader->get_sequence_count(),
			params->filter_k, 
			&matrix );
	std::cerr << "ok" << std::endl;

	std::cerr << "   Index-table for kDP alignments...";
	_4mer_table = new Recycle_table( params->max_number_of_4mer_table_representants, 
			params->kdp_k, 
			&matrix );
	std::cerr << "ok" << std::endl;

	if (params->kmer_matrices){
		std::cerr << "   4-mer matrix...";
		_4mer_matrix = new Kmer_matrix(4, matrix.get_scale());
		//	    _4mer_matrix->fill_kmer_matrix("/archive/user/maria/data_backup/kmer_matrices/matrix_4mer_db_nr30_2008_neff2.5_hhfilter_id50_maxseqnum_50_new_indexes.dat");
		//	    _4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_db_nr30_2008_neff2.5_hhfilter_id50_kalign_alns_mincount5_maxseqnum_80.dat");
		//	    _4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_XX0XX_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_50.dat");
		//_4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_50.dat");
		_4mer_matrix->fill_kmer_matrix(params->_4mer_m_file);
		std::cerr << "ok" << std::endl;

		std::cerr << "   2-mer matrix...";
		_2mer_matrix = new Kmer_matrix(2, matrix.get_scale());
		//	    _2mer_matrix->fill_kmer_matrix("/archive/user/maria/data_backup/kmer_matrices/matrix_2mer_db_nr30_2008_neff2.5_hhfilter_id50_maxseqnum_50_new_indexes.dat");
		// _2mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_2mer_db_nr30_2008_neff2.5_hhfilter_id50_kalign_alns_maxseqnum_50.dat");
		//		 _2mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_2mer_X0X_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_10.dat");
		//_2mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_2mer_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_50.dat");
		_2mer_matrix->fill_kmer_matrix(params->_2mer_m_file);
		std::cerr << "ok" << std::endl;

	}

	std::cerr << "   Filter k-mers...";
	_6mer       = new Kmer( params->filter_k, 
			params->filter_kmer_similarity_threshold,
			params->score_correction,
			_4mer_matrix,
			_2mer_matrix);
	std::cerr << "ok (k=" << params->filter_k << ", thr=" << params->filter_kmer_similarity_threshold << ")" << std::endl;

	std::cerr << "   kDP k-mers...";
	_4mer       = new Kmer( params->kdp_k, 
			params->kdp_kmer_similarity_threshold,
			params->score_correction,
			_4mer_matrix);
	std::cerr << "ok (k=" << params->kdp_k << ", thr=" << params->kdp_kmer_similarity_threshold << ")" <<  std::endl;

	_3mer       = new Kmer( 3,
			params->filter_kmer_similarity_threshold,
			params->score_correction,
			_4mer_matrix,
			_2mer_matrix);

	std::cerr << "   kDP...";
	if( params->check_coverage_seq_len )
		_4meraln    = new Kmeraln( params->kdp_k, 
				&matrix, params->kdp_p_m, 
				params->kdp_delta, 
				params->kdp_G, 
				params->kdp_E, 
				params->kdp_F, 
				params->kdp_Lamda, 
				params->kdp_K, 
				params->kdp_H,
				db_reader->get_character_count(),
				db_reader->get_sum_log_template_len(),
				(float)db_reader->get_sequence_count(),
				2,
				params->coverage_seq_len-0.05 );
	else
		_4meraln    = new Kmeraln( params->kdp_k, 
				&matrix, 
				params->kdp_p_m, 
				params->kdp_delta, 
				params->kdp_G, 
				params->kdp_E, 
				params->kdp_F,
				params->kdp_Lamda, 
				params->kdp_K, 
				params->kdp_H,
				db_reader->get_character_count(),
				db_reader->get_sum_log_template_len(),
				(float)db_reader->get_sequence_count() );
	std::cerr << "ok" << std::endl;

	offset_memory_usage = get_memory_usage();
	//std::cerr << "   Current memory         : " << offset_memory_usage/(1024*1024) << "MB" << std::endl;
	//std::cerr << "   Available memory       : " << params->memory_limit/(1024*1024) << "MB" << std::endl;

/*	std::cout << "offset memory usage: " << offset_memory_usage << "\n";
	std::cout << "+ " << (1024*1024) << "\n";
	std::cout << "= " << (offset_memory_usage+1024*1024) << "\n";
	std::cout << "memory limit: " << params->memory_limit <<"\n";*/

	if( (offset_memory_usage+1024*1024)>params->memory_limit ){ //if( (offset_memory_usage+50*1024*1024)>params->memory_limit ){
		int min_mem = (offset_memory_usage/(1024*1024)) + 50;
		throw MyException("Memory limit too low! >%i MB required.", min_mem);
	}

	double seq_mem_limit = (double)params->memory_limit - (double)offset_memory_usage;
//	std::cout << "seq_mem_limit: " << seq_mem_limit << "\n";

	double estimated_mem = (sizeof(Base_table::_simple_rkmer)+1)*db_reader->get_character_count();
	estimated_mem       += db_reader->get_sequence_count()*5*(sizeof(size_t));
	if (params->profile_query) estimated_mem += db_reader->get_sequence_count() * 20 * sizeof(float);
//	std::cout << "estimated_mem: " << estimated_mem << "\n";

	if( estimated_mem<seq_mem_limit ){
		std::cerr << "   No swapping required." << std::endl;
		selfcompare_mem_limit = params->memory_limit;
	}else{
		std::cerr << "   Memory limit will be reached -> swapping." << std::endl;
		selfcompare_mem_limit = (size_t)(0.7*seq_mem_limit + offset_memory_usage);
	}
	block_index          = 1;
}

size_t K_Clust::get_memory_usage(){
	size_t mem = clusters->get_memory_usage();
	mem += block_man->get_memory_usage();
	mem += _6mer->get_memory_usage();
	mem += _4mer->get_memory_usage();
	mem += _4meraln->get_memory_usage();
	mem += db_reader->get_memory_usage();
	mem += table->get_memory_usage();
	mem += _4mer_table->get_memory_usage();
	mem += list_of_new_candidates.size() * 2 *sizeof(void*);
	return mem;
}

void K_Clust::prepare(const size_t width) throw (std::exception){
	string aln_dir = params->get_aln_dir_old();
	//	string hhconslogf = params->get_input_dir() + "hhconsensus.log";
	//	string hhmakelogf = params->get_input_dir() + "hhmake.log";

	int ret = 0;
	// copy old db_sorted file in the new directory
	ret = system(("cp " + params->get_dbsortedfile_old() + " " + params->get_working_dir()).c_str());
	if (ret != 0) throw MyException ("Cannot execute '%s'", ("cp " + params->get_dbsortedfile_old() + " " + params->get_working_dir()).c_str());

	DIR *dir = opendir(aln_dir.c_str());
	if (dir == NULL) throw MyException ("Cannot open the directory with multiple alignments '%s'", aln_dir.c_str());

	struct stat stFileInfo;
	struct dirent * entry;

	int idx = 0;

	std::ofstream out ((params->get_dbfile()).c_str());
	std::ifstream in_seq;

	char* hbuffer = new char [READ_HBUFFER_SIZE];
	char* sbuffer = new char [READ_SBUFFER_SIZE];

	while((entry = readdir(dir)) != NULL)
	{
		if (string(entry->d_name) == "." || string(entry->d_name) == "..") continue;

		std::cerr << string(entry->d_name) << "  ";
		string dir1_path = aln_dir + string(entry->d_name) + "/";
		DIR *dir1 = opendir(dir1_path.c_str());

		struct dirent * entry1;
		while ((entry1 = readdir(dir1)) != NULL){
			string fname = string(entry1->d_name);

			if (fname == "." || fname == "..") continue;
			if(fname.substr(fname.length()-4, fname.length()).compare(".a3m") != 0) continue;

			// write consensus sequence with representative index as the header into the temporary
			// database file
			fname = fname.substr(0, fname.length()-4);
			idx = Converter::get_idx_for_file(fname.c_str());
			out << ">" << idx << "_consensus\n";

			in_seq.open((dir1_path + fname + ".a3m").c_str());
			in_seq.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
			if (hbuffer[0] == '#') in_seq.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
			in_seq.getline(sbuffer, READ_SBUFFER_SIZE, '\n');

			//write the sequence into a temporary database file
			bool endofseq = false;
			for (size_t i = 0; (i<READ_SBUFFER_SIZE) && (!endofseq); i+=width){
				if (sbuffer[i]=='\0')
					break;
				for (size_t j = i; j<(i+width); j++){
					if (sbuffer[j]=='\0'){
						endofseq = true;
						break;
					}
					out << sbuffer[j];
				}
				out << '\n';
			}
			in_seq.close();

		}
		closedir(dir1);

	}
	std::cerr << "\n";
	closedir(dir);

	out.clear();
	out.close();

	delete[] hbuffer;
	delete[] sbuffer;

}

void K_Clust::pref_benchmark_scop_tpfp() throw (std::exception){
	_init();
	table->reset();

	int kmer_count;
	float _6mer_per_pos;
	double kmer6_av = 0.0;
	int seqs_6 = 0;

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);
		table->add_representative_spaced( s );
		block_man->save_representative( s );
	}

	db_reader->reset();
	int num_matches = 0;
	int num_seqs = 0;

	float max_score = 0.0f;

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);

		std::string q_id (s->get_header());
		q_id = q_id.substr(1, 7);
		std::ofstream outfile ((params->get_working_dir() + q_id).c_str());
		if( !outfile ) throw MyException("Cannot open '%s' for writing!", (params->get_working_dir() + q_id).c_str());

		if (params->count_idents){
			kmer_count  = _6mer->create_idents_list(s);
			table->match( _6mer );
		}
		else {
			kmer_count  = _6mer->create_kmer_list_fast_spaced(s);
			table->match( _6mer );
		}

		if (s->length() >= 6){
			_6mer_per_pos = ((float) kmer_count) / ((float) (s->length() - 5));
			kmer6_av += _6mer_per_pos;
			seqs_6++;
		}
		// sequence length, score threshold, coverage
		Base_table::MatchIterator it  = table->begin( (float)s->length(), 0.0, 0.0 );
		Base_table::MatchIterator end = table->end();

		// there are some representative sequences above filtering threshold
		while( it!=table->end()){
			Sequence *r = block_man->get_representative( it.get_sequence_index() );
			if (s->index() == r->index()){
				++it;
				continue;
			}
			num_matches++;
			std::string t_id (r->get_header());
			t_id = t_id.substr(1, 7);
			if (max_score < it.get_score()) max_score = it.get_score();
			outfile << t_id << "\t" << it.get_score() << "\n";
			++it;

		}
		outfile.close();
		num_seqs++;
		if (num_seqs % 1000 == 0) std::cout << num_seqs << " query sequences processed.\n";
	}
	std::cout << "matches: " << num_matches << "\n";
	std::cout << "max_score: " << max_score << "\n";
	std::cout << kmer6_av << " " << seqs_6 << "\n";
	std::cout << "6-mers per pos: " << kmer6_av/(double)seqs_6 << "\n";
	/*	string rm1 = "rm -r " +  params->get_working_dir() + "tmp";
		string rm2 = "rm " +  params->get_working_dir() + "db_sorted.fas";
		std::cout << rm1 << "\n";
		std::cout << rm2 << "\n";
		system(rm1.c_str());
		system(rm2.c_str());*/
}

void K_Clust::pref_benchmark_scop_scatterplot() throw (std::exception){
	_init();
	table->reset();

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);
		table->add_representative_spaced( s );
		block_man->save_representative( s );
	}

	db_reader->reset();
	int num_matches = 0;
	int num_seqs = 0;

	float max_score = 0.0f;
	char header [1000];

	std::string q_id;
	char q_class;
	int q_fold;
	int q_sf;

	std::string t_id;
	char t_class;
	int t_fold;
	int t_sf;

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);

		strcpy(header, s->get_header());

		// scop id
		q_id = strtok(header," ");
		q_id = q_id.substr(1, 7);
		// class
		q_class = strtok(NULL,".")[0];
		// fold
		q_fold = atoi(strtok(NULL,"."));
		// superfamily
		q_sf = atoi(strtok(NULL,"."));
		std::ofstream outfile ((params->get_working_dir() + q_id).c_str());
		if( !outfile ) throw MyException("Cannot open '%s' for writing!", (params->get_working_dir() + q_id).c_str());

		if (params->count_idents){	
			_3mer->create_idents_list( s );
			table->match( _3mer );
		}
		else{
			_6mer->create_kmer_list_fast_spaced( s );
			table->match( _6mer );
		}

		// sequence length, score threshold, coverage
		Base_table::MatchIterator it  = table->begin( (float)s->length(), 0.0, 0.0 );
		Base_table::MatchIterator end = table->end();

		// there are some representative sequences above filtering threshold
		while( it!=table->end()){
			Sequence *r = block_man->get_representative( it.get_sequence_index() );
			if (s->index() == r->index()){
				++it;
				continue;
			}
			strcpy(header, r->get_header());

			// scop id
			t_id = strtok(header," ");
			t_id = t_id.substr(1, 7);
			// class
			t_class = strtok(NULL,".")[0];
			// fold
			t_fold = atoi(strtok(NULL,"."));
			// superfamily
			t_sf = atoi(strtok(NULL,"."));

			if (q_class == t_class && q_fold == t_fold && q_sf == t_sf){

				if (max_score < it.get_score()/(float)s->length()) max_score = it.get_score()/(float)s->length();
				outfile << t_id << "\t" << it.get_score()/(float)s->length() << "\n";
				num_matches++;
			}
			++it;

		}
		outfile.close();
		num_seqs++;
		if (num_seqs % 1000 == 0) std::cout << num_seqs << " query sequences processed.\n";
	}
	/*	string rm1 = "rm -r " +  params->get_working_dir() + "tmp";
		string rm2 = "rm " +  params->get_working_dir() + "db_sorted.fas";
		std::cout << rm1 << "\n";
		std::cout << rm2 << "\n";
		system(rm1.c_str());
		system(rm2.c_str());*/
	std::cout << "matches: " << num_matches << "\n";
	std::cout << "max_score: " << max_score << "\n";
}

void K_Clust::kDP_benchmark_scop_tpfp() throw (std::exception){
	_init();
	_4mer_table->reset();

	int kmer_count;
	float _4mer_per_pos;
	float kmer4_av = 0.0f;
	int seqs_4 = 0;

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);
		_4mer_table->add_representative_spaced( s );
		block_man->save_representative( s );
	}

	db_reader->reset();
	int num_matches = 0;
	int num_seqs = 0;

	float max_score = 0.0f;

	while(db_reader->has_next() ){
		Sequence *s = db_reader->get_next(false, false);

		std::string q_id (s->get_header());
		q_id = q_id.substr(1, 7);
		std::ofstream outfile ((params->get_working_dir() + q_id).c_str());

		if( !outfile ) throw MyException("Cannot open '%s' for writing!", (params->get_working_dir() + q_id).c_str());

		kmer_count = _4mer->create_kmer_list_fast_spaced(s);
		if (s->length() >= 4){
			_4mer_per_pos = ((float) kmer_count) / ((float) (s->length() - 3));
			kmer4_av += _4mer_per_pos;
			seqs_4++;
		}
		_4mer_table->match_full( _4mer );

		Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
		Recycle_table::MatchIterator alns_end = _4mer_table->end();
		float score=0.0f;

		while( alns_it!=alns_end ){
			Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );
			if (s->index() == r->index()){
				++alns_it;
				continue;
			}
			std::string t_id (r->get_header());
			t_id = t_id.substr(1, 7);

			num_matches++;

			_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r );

			outfile << t_id << "\t" << _4meraln->get_score() << "\n";
			if (_4meraln->get_score() > max_score) max_score = _4meraln->get_score();

			++alns_it;
		}

		outfile.close();
		num_seqs++;
		if (num_seqs % 1000 == 0) std::cout << num_seqs << " query sequences processed.\n";
	}
	std::cout << "matches: " << num_matches << "\n";
	std::cout << "max_score: " << max_score << "\n";
	std::cout << "4-mers per pos: " << kmer4_av/(float)seqs_4 << "\n";
	/*	string rm1 = "rm -r " +  params->get_working_dir() + "tmp";
		string rm2 = "rm " +  params->get_working_dir() + "db_sorted.fas";
		std::cout << rm1 << "\n";
		std::cout << rm2 << "\n";
		system(rm1.c_str());
		system(rm2.c_str());*/
}

void K_Clust::kDP_benchmark_scop_TMscore() throw (std::exception){
	if (params->kmer_matrices){
		std::cerr << "   4-mer matrix...";
		_4mer_matrix = new Kmer_matrix(4, matrix.get_scale());
		//	    _4mer_matrix->fill_kmer_matrix("/archive/user/maria/data_backup/kmer_matrices/matrix_4mer_db_nr30_2008_neff2.5_hhfilter_id50_maxseqnum_50_new_indexes.dat");
		//		_4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_db_nr30_2008_neff2.5_hhfilter_id50_kalign_alns_mincount5_maxseqnum_80.dat");
		_4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_XX0XX_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_50.dat");
		//	    _4mer_matrix->fill_kmer_matrix("/cluster/user/maria/kClust_data/kmer_matrices/matrix_4mer_nr20f_151208_neff2.5_hhblits_last_round_maxseqnum_50.dat");
		std::cerr << "ok" << std::endl;
	}

	std:string input_dir = params->get_input_dir();
	params->profile_query = false;

	DIR *dir = opendir(input_dir.c_str());
	if (dir == NULL) throw MyException ("Cannot open directory '%s'", input_dir.c_str());

	struct stat stFileInfo;
	struct dirent * entry;

	int count = 0;

	while((entry = readdir(dir)) != NULL)
	{

		string fname = string(entry->d_name);

		if (fname.length() < 4 ) continue;
		if(fname.substr(fname.length()-3, fname.length()).compare("_db") != 0) continue;

		if (count > 0){
			delete db_reader;
			delete clusters;
			delete block_man;
			delete table;
			delete _4mer_table;
			delete _4mer;
			delete _4meraln;
		}

		std::cerr << "======================================================" << std::endl;
		params->set_dbfile(input_dir + entry->d_name);
		std::cout << "Current database file : " << input_dir << entry->d_name << "\n";

		Fasta_db_reader dbr( params->get_dbfile(), &matrix, false, false );
		dbr.seek_sort( params->get_dbsortedfile() );
		std::cout << "# sequences: " << dbr.get_sequence_count() << "\n";

		db_reader   = new Fasta_db_reader( params->get_dbsortedfile(), &matrix, params->score_correction, params->profile_query, params->get_aln_dir_old(), params->get_clustersfile_old(), params->get_headersfile_old() );
		clusters    = new Clusters( db_reader->get_sequence_count(), &matrix );
		block_man   = new BlockManager( db_reader, clusters );
		table       = new Base_table( db_reader->get_sequence_count(),
				params->filter_k,
				&matrix );
		_4mer_table = new Recycle_table( params->max_number_of_4mer_table_representants,
				params->kdp_k,
				&matrix );
		_4mer       = new Kmer( params->kdp_k,
				params->kdp_kmer_similarity_threshold,
				params->score_correction,
				_4mer_matrix);

		if( params->check_coverage_seq_len )
			_4meraln    = new Kmeraln( params->kdp_k,
					&matrix, params->kdp_p_m,
					params->kdp_delta,
					params->kdp_G,
					params->kdp_E,
					params->kdp_F,
					params->kdp_Lamda,
					params->kdp_K,
					params->kdp_H,
					db_reader->get_character_count(),
					db_reader->get_sum_log_template_len(),
					(float)db_reader->get_sequence_count(),
					1,
					params->coverage_seq_len-0.05 );
		else
			_4meraln    = new Kmeraln( params->kdp_k,
					&matrix,
					params->kdp_p_m,
					params->kdp_delta,
					params->kdp_G,
					params->kdp_E,
					params->kdp_F,
					params->kdp_Lamda,
					params->kdp_K,
					params->kdp_H,
					db_reader->get_character_count(),
					db_reader->get_sum_log_template_len(),
					(float)db_reader->get_sequence_count() );


		_4mer_table->reset();

		while(db_reader->has_next() ){
			Sequence *s = db_reader->get_next(false, false);
			_4mer_table->add_representative_spaced( s );
			block_man->save_representative( s );
		}

		db_reader->reset();
		int num_matches = 0;
		int num_seqs = 0;

		float max_score = 0.0f;

		char header [1000];

		std::string q_id;
		char q_class;
		int q_fold;
		int q_sf;

		std::string t_id;
		char t_class;
		int t_fold;
		int t_sf;

		while(db_reader->has_next() ){
			Sequence *s = db_reader->get_next(false, false);

			strcpy(header, s->get_header());

			// scop id
			q_id = strtok(header," ");
			q_id = q_id.substr(1, 7);
			// class
			q_class = strtok(NULL,".")[0];
			// fold
			q_fold = atoi(strtok(NULL,"."));
			// superfamily
			q_sf = atoi(strtok(NULL,"."));

			std::ofstream outfile ((params->get_working_dir() + q_id + ".kdp").c_str());

			if( !outfile ) throw MyException("Cannot open '%s' for writing!", (params->get_working_dir() + q_id + ".kdp").c_str());

			_4mer->create_kmer_list_fast_spaced(s);
			_4mer_table->match_full( _4mer );

			Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
			Recycle_table::MatchIterator alns_end = _4mer_table->end();
			float score=0.0f;

			while( alns_it!=alns_end ){
				Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );

				num_matches++;

				_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r );

				if (_4meraln->get_score() > max_score) max_score = _4meraln->get_score();

				_4meraln->print_alignment(outfile);
				outfile << "\n\n";
				//				}

				++alns_it;
			}

			outfile.close();
			num_seqs++;
			if (num_seqs % 1000 == 0) std::cout << num_seqs << " query sequences processed.\n";
		}

		count++;
		std::cout << count << " families processed\n";
	}
}


void K_Clust::cluster() throw (std::exception){
	std::cerr << "======================================================" << std::endl;
	counter__all_kmer_alns         = 0;
	counter__kmeraln_trigger_seqs  = 0;
	counter__max_kmeraln_count     = 0;
	_init();

	std::cerr << "Clustering..." << std::endl;
	std::cerr << "Clustering threshold: " << params->clustering_threshold << "\n";

	BlockManager::Block b( params->get_tmp_dir().c_str(), block_index++ );

	_selfcompare_db_reader( b );
	block_man->add_block( b );

	b.print(std::cerr, "   >");
	// write index table to file, representative sequences remain in the memory in block manager
	if( params->refinement || db_reader->has_next() ){
		std::cerr << "   Swapping index-table ('" << b.get_filename() << "') ...";
		table->write_table_to_file( b.get_filename().c_str() );
		std::cerr << "ok" << std::endl;
	}
	while( db_reader->has_next() ){
		std::cerr << " Next column ..." << std::endl;
		b.set_index( block_index++ );

		BlockManager::Block skip_block = block_man->get_current_table_block();
		std::cerr << "  Reading sequences and comparing to segment ";
		std::cerr << block_man->get_current_table_block().get_index() << std::endl;

		_compare_db_reader( b );

		b.print( std::cerr, "   >");

		BlockManager::Iterator it = block_man->begin();
		while( it!=block_man->end() ){

			if( *it==skip_block ){ 
				//std::cerr << " Skipping segment " << (*it).get_index() <<  std::endl;
				++it;
				continue;
			}

			const float minimum_coverage = b.length_of_longest_seq/(float)(*it).length_of_shortest_seq;
			if( minimum_coverage<params->coverage_seq_len ){
				std::cerr << "  Segment " << (*it).get_index() << " contains no sequences in coverage range";
				std::cerr << " (min:" <<(*it).length_of_shortest_seq << "|max:" << (*it).length_of_longest_seq << ")" << std::endl;
				++it;
				continue;
			}

			std::cerr << "  Comparing sequences to segment "<< (*it).get_index() << std::endl;
			std::cerr << "   Loading '"<< (*it).get_filename() << "' ... ";
			table->load_from_file( (*it).get_filename().c_str() );
			std::cerr << "ok" << std::endl;

			std::cerr << "   Loading representatives ...";
			block_man->load( (*it) );
			std::cerr << "ok" << std::endl;

			_compare_list();

			++it;	
		}
		_selfcompare_list( b );
		block_man->add_block( b );
		b.print( std::cerr, "   >" );

		std::cerr << "   Swapping index-table ('" << b.get_filename() << "') ...";
		table->write_table_to_file( b.get_filename().c_str() );
		std::cerr << "ok" << std::endl;

		if( block_man->get_current_table_block().end!=db_reader->get_sequence_count()+1 ){
			db_reader->move_to( b.end );
		}else{
			std::cerr << "All sequences processed!" << std::endl;
			break;
		}
	}

	std::cerr << "======================================================"                                                 << std::endl;
	std::cerr << "Number of clusters                  : " << clusters->get_number_of_clusters()                           << std::endl;
	std::cerr << "======================================================" << std::endl;
	std::cout << "Average 6mer number per position: " << kmer6_av/(float)seqs_6 << "\n";
	std::cout << "Average 4mer number per position: " << kmer4_av/(float)seqs_4 << "\n";
	std::cout << "Average 6mer db matches per sequence: " << matches_6mer/(float)seqs_6 << "\n";
	std::cout << "Average 4mer db matches per sequence: " << matches_4mer/(float)seqs_6 << "\n";
	std::cout << "# Sequences in kDP: " << counter__kmeraln_trigger_seqs << "\n";
	std::cout << "# of kmer alns: " << counter__all_kmer_alns << "\n";

	if (params->write_time_bench){
		std::ofstream out_pref ((params->get_working_dir() + "tb_prefiltering.dat").c_str());
		std::ofstream out_kdp ((params->get_working_dir() + "tb_kDP.dat").c_str());
		std::ofstream out_seq ((params->get_working_dir() + "tb_all.dat").c_str());

		out_pref << "# seq_ID  seq_length  6mers_per_position  time  time/seq_length\n";
		out_kdp << "# seq_ID  seq_length  4mers_per_position  alignments  time  time/(seq_length^2)\n";
		out_seq << "# seq_ID  seq_length  6mers_per_position  4mers_per_position  alignments  time_prefiltering+time_kDP  time_prefiltering/seq_length+time_kDP/(seq_length^2)\n";


		while (q_pref->size() > 0){
			time_benchmark* t = q_pref->top();
			q_pref->pop();
			out_pref << t->seq_id << "  " << t->length << "  " << t->_6mer_per_pos << "  " << t->time << "  " << t->time_per_pos << "\n";
			Sequence* s = db_reader->get(t->seq_id, true, false, true);
			out_pref << s->get_header() << "\n";
			out_pref << s->sequence() << "\n";
			delete t;
			delete s;
		}
		out_pref.close();

		while (q_kdp->size() > 0){
			time_benchmark* t = q_kdp->top();
			q_kdp->pop();
			out_kdp << t->seq_id << "  " << t->length << "  " << t->_4mer_per_pos << "  " << t->aln_num << "  " << t->time << "  " << t->time_per_pos << "\n";
			Sequence*  s = db_reader->get(t->seq_id, true, false, true);
			out_kdp << s->get_header() << "\n";
			out_kdp << s->sequence() << "\n";
			delete t;
			delete s;
		}
		out_kdp.close();

		while (q_seq->size() > 0){
			time_benchmark* t = q_seq->top();
			q_seq->pop();
			out_seq << t->seq_id << "  " << t->length << "  " << t->_6mer_per_pos << "  " << t->_4mer_per_pos << "  " << t->aln_num << "  " << t->time << "  " << t->time_per_pos << "\n";
			Sequence* s = db_reader->get(t->seq_id, true, false, true);
			out_seq << s->get_header() << "\n";
			out_seq << s->sequence() << "\n";
			delete t;
			delete s;
		}
		out_seq.close();
	}
}

void K_Clust::refine() throw (std::exception){

	std::cerr << "======================================================" << std::endl;

	counter__all_kmer_alns         = 0;
	counter__kmeraln_trigger_seqs  = 0;
	counter__max_kmeraln_count     = 0;


	std::cerr << "Removing members...";
	clusters->remove_members();
	std::cerr << "ok" << std::endl;;

	std::cerr << "Starting refinement..." << std::endl; 
	BlockManager::Iterator b_it = block_man->begin();
	while( b_it!=block_man->end() ){
		BlockManager::Block current_block = *b_it;
		std::cerr << " Loading table '"<< current_block.get_filename() << "' ...";
		table->load_from_file( current_block.get_filename().c_str() );
		std::cerr << "ok" << std::endl;
		std::cerr << " Loading corresponding representative sequences ...";
		block_man->load( current_block );
		std::cerr << "ok" << std::endl;
		std::cerr << " Comparing member sequences to table " << current_block.get_index() << std::endl;
		db_reader->reset();
		size_t mem_seq_counter=0;
		while( db_reader->has_next() ){

			if( clusters->is_representative(db_reader->get_index_of_next_sequence())) continue;
			Sequence *s = db_reader->get_next(false, false);

			//check if sequences are too short -> there will be no more relevant sequences for tis block
			const float minimum_cov = s->length()/(float)current_block.length_of_shortest_seq;
			if( minimum_cov<params->coverage_seq_len ){
				delete s;
				break;
			}

			//check if the sequence is too long -> go on until sequences with length  in coverage range are processed
			const float current_cov = current_block.length_of_longest_seq/(float)s->length();
			if( current_cov<params->coverage_seq_len ){
				delete s;
				continue;
			} 

			if( ++mem_seq_counter%params->sequence_block_size==0 ) std::cerr << ".";

			if (! params->profile_query) _6mer->create_kmer_list_fast_spaced( s );
			else _6mer->create_kmer_list_spaced(s);
			table->match( _6mer );

			Base_table::MatchIterator it  = table->begin(   (float)s->length(),
					params->kmer_score_threshold,
					params->coverage_seq_len );
			Base_table::MatchIterator end = table->end();
			if( it!=end ){
				++counter__kmeraln_trigger_seqs; 
				_4mer_table->reset();
				if (! params->profile_query) _4mer->create_kmer_list_fast_spaced( s );
				else _4mer->create_kmer_list_spaced(s);

				size_t count_homologs=0;
				while( it!=end && count_homologs<params->max_number_of_4mer_table_representants ){
					++counter__all_kmer_alns;
					Sequence *r = block_man->get_representative( it.get_sequence_index() );
					_4mer_table->add_representative_spaced( r );
					++it;
					++count_homologs;
				}

				if( count_homologs>counter__max_kmeraln_count ){
					counter__max_kmeraln_count = count_homologs;
				}

				if( it!=table->end() && count_homologs==params->max_number_of_4mer_table_representants )
					throw MyException("Sequence limit of 4-mer table reached! (%zu) ", params->max_number_of_4mer_table_representants);


				_4mer_table->match_full( _4mer );
				Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
				Recycle_table::MatchIterator alns_end = _4mer_table->end();
				size_t idx  = 0;
				float score = 0.0f;
				// for each sequence aligned to the query sequence there is a Recycle_table of matches. MatchIterators iterate over these Recycle_tables.
				while( alns_it!=alns_end ){
					Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );

					_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r);
					float current_score = _4meraln->get_rel_score();
					//					if( current_score>score && _4meraln->get_aln_cov_short()>=params->aln_coverage_short && _4meraln->get_aln_cov_long()>=params->aln_coverage_long ){
					if( current_score>score && _4meraln->get_aln_cov_long()>=params->aln_coverage_long ){
						score = current_score;
						idx = alns_it.get_sequence_index();
					}
					++alns_it;
				}
				if( idx>0 ){
					clusters->add_member( idx, s->index(), score );
				}
			}else{
				//The member has no match with any representative in this table
			}
			delete s;
		}
		std::cerr << std::endl;
		++b_it;
	}


	std::cerr << "======================================================"                                                 << std::endl;
	std::cerr << "Reassigning member-sequences...";
	clusters->evaluate_members();
	std::cerr << "ok" << std::endl;

	std::cerr << "======================================================"                                                 << std::endl;
	std::cerr << "Number of clusters                  : " << clusters->get_number_of_clusters()                           << std::endl;
	std::cerr << "======================================================" << std::endl;
}



void K_Clust::_selfcompare_db_reader( BlockManager::Block &b ) throw (std::exception){
	std::cout << "Function: _selfcompare_db_reader\n";
	size_t mem       = 0;
	table->reset();
	b.begin          = db_reader->get_index_of_next_sequence();
	std::cout << "begin of the block: " <<b.begin <<"\n";
	size_t minlength = INT_MAX;
	size_t maxlength = 0;
	std::cerr << "   ";
	while( mem<selfcompare_mem_limit && db_reader->has_next() ){
		// dots
		std::cerr << ".";
		size_t c=0;
		while( c < params->sequence_block_size && db_reader->has_next() ){
			Sequence *s = db_reader->get_next(false, false);
			++c;
			__selfcompare( s, minlength, maxlength, b.get_index());
		}
		mem = get_memory_usage();
	}
	//new line behind the dots
	std::cerr << std::endl;
	b.end                    = db_reader->get_index_of_next_sequence();
	b.length_of_shortest_seq = minlength;
	b.length_of_longest_seq  = maxlength;
	b.characters             = table->get_number_of_rkmers();
}

void K_Clust::_selfcompare_list( BlockManager::Block &b ) throw (std::exception){
	std::cout << "Function: _selfcompare_list\n";
	table->reset();
	std::list<Sequence*>::iterator it  = list_of_new_candidates.begin();
	std::list<Sequence*>::iterator end = list_of_new_candidates.end();
	size_t minlength                   = INT_MAX;
	size_t maxlength                   = 0;
	size_t c                           = 0;
	std::cerr << "   ";
	while( it!=end ){
		Sequence *s = *it;
		__selfcompare( s , minlength, maxlength, b.get_index() );
		++it;
		++c;
		if( c%params->sequence_block_size==0 ) std::cerr << ".";
	}
	//new line behind the dots
	std::cerr << std::endl;
	list_of_new_candidates.clear();
	b.length_of_shortest_seq = minlength;
	b.length_of_longest_seq  = maxlength;
	b.characters = table->get_number_of_rkmers();
}

void K_Clust::__selfcompare( Sequence *s, size_t &minlength, size_t &maxlength, const size_t table_idx) throw (std::exception){

#ifdef __CM_DEBUG_PRINT_LETTERS
	std::cerr << "l" << s->length() << "\n";
#endif

	size_t kmer_count = 0;
	//	std::cout << "\n";
	//	std::cout << (db_reader->get_idx2orig())[s->index()] << "\n";
	//	std::cout << s->get_header() << "\n";
	//	std::cout << "Selfcomparing sequence " << s->index() << "\n";

	Chronometer chron;
	double ms_prefilter = 0.0;
	double ms_kdp = 0.0;
	float _6mer_per_pos = 0.0;
	float _4mer_per_pos = 0.0;
	chron.start();

	if (! params->profile_query) kmer_count  = _6mer->create_kmer_list_fast_spaced( s );
	else kmer_count = _6mer->create_kmer_list_spaced(s);
	if (s->length() >= 6){
		_6mer_per_pos = ((float) kmer_count) / ((float) (s->length() - 5));
		kmer6_av += _6mer_per_pos;
		seqs_6++;
	}


	matches_6mer += table->match( _6mer );

	chron.stop();
	ms_prefilter = 1000.0*chron.getSeconds();

	time_benchmark* t = new time_benchmark(s->index(), s->length(), _6mer_per_pos, 0.0f, 0,  ms_prefilter, ms_prefilter/s->length());
	insert_into_queue(q_pref, t);

	int kmer_alns = 1;
	chron.start();
	Base_table::MatchIterator it  = table->begin( 	(float)s->length(),
			params->kmer_score_threshold,
			params->coverage_seq_len );
	Base_table::MatchIterator end = table->end();

	// there are some representative sequences above filtering threshold
	int alns = 0;
	if( it!=end ){
		++counter__kmeraln_trigger_seqs;
		_4mer_table->reset();

		if (! params->profile_query) kmer_count  = _4mer->create_kmer_list_fast_spaced( s );
		else kmer_count  = _4mer->create_kmer_list_spaced(s);
		if (s->length() >= 4){
			_4mer_per_pos = ((float) kmer_count) / ((float) (s->length() - 3));
			kmer4_av += _4mer_per_pos;
			seqs_4++;
		}

		// add representatives from this block which passed the kmer filtering step to the Recycle_table
		// Recycle_table stores position i in seq. x and j in  seq. y for each match
		size_t count_homologs=0;
		while( it!=table->end() && count_homologs < params->max_number_of_4mer_table_representants){
			++counter__all_kmer_alns;
			Sequence *r = block_man->get_representative( it.get_sequence_index() );
			_4mer_table->add_representative_spaced( r );
			++it;
			++count_homologs;
		}

		if( count_homologs>counter__max_kmeraln_count ){
			counter__max_kmeraln_count = count_homologs;
		}

		//		if( it!=table->end() && count_homologs==params->max_number_of_4mer_table_representants )
		//			throw MyException("Sequence limit of 4-mer table reached! (%i) ", params->max_number_of_4mer_table_representants);

		matches_4mer += _4mer_table->match_full( _4mer );
		Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
		Recycle_table::MatchIterator alns_end = _4mer_table->end();
		size_t idx=0;
		float score=0.0f;


		// compute alignments
		while( alns_it!=alns_end ){
			alns++;
			Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );
			//			std::cout << "Aligning " << s->index() << " and " << r->index() << "\n";
			_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r );
			//			std::cout << "\n";

			//			_4meraln->print_alignment(std::cout);

			float current_score = _4meraln->get_rel_score();
			float current_evalue = _4meraln->get_evalue();

			//			if(  current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_short()>=params->aln_coverage_short && _4meraln->get_aln_cov_long()>=params->aln_coverage_long){
			if(  current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_long()>=params->aln_coverage_long){
				score = current_score;
				idx = alns_it.get_sequence_index();
				break;
			}
			++alns_it;
		}


		// sequence s is not a representative, add to cluster
		if( idx>0 ){
			clusters->add_member( idx, s->index(), score );

			chron.stop();
			ms_kdp = 1000.0*chron.getSeconds();
			t = new time_benchmark(s->index(), s->length(), 0.0f, _4mer_per_pos, alns, ms_kdp, ms_kdp/(s->length()*s->length()));
			insert_into_queue(q_kdp, t);

			t = new time_benchmark(s->index(), s->length(), _6mer_per_pos, _4mer_per_pos, alns, ms_prefilter + ms_kdp, ms_prefilter/s->length() + ms_kdp/(s->length()*s->length()));
			insert_into_queue(q_seq, t);

			delete s;
			return;
		}

		chron.stop();
		ms_kdp= 1000.0*chron.getSeconds();

		t = new time_benchmark(s->index(), s->length(), 0.0f, _4mer_per_pos, alns, ms_kdp, ms_kdp/(s->length()*s->length()));
		insert_into_queue(q_kdp, t);
	}

	t = new time_benchmark(s->index(), s->length(), _6mer_per_pos, _4mer_per_pos, alns, ms_prefilter + ms_kdp, ms_prefilter/s->length() + ms_kdp/(s->length()*s->length()));
	insert_into_queue(q_seq, t);

	// sequence s is a representative
	clusters->create_cluster_for( s->index(), table_idx);
	table->add_representative_spaced( s );
	block_man->save_representative( s );
	minlength = std::min( s->length(), minlength );
	maxlength = std::max( s->length(), maxlength );
}

void K_Clust::_compare_db_reader( BlockManager::Block &b ){
	std::cout << "Function: compare_db_reader\n";
	size_t mem               = 0;
	size_t seq_mem           = 0;
	size_t chars             = 0;
	const float shortest     = block_man->get_current_table_block().length_of_shortest_seq;
	b.begin                  = db_reader->get_index_of_next_sequence();
	b.length_of_shortest_seq = INT_MAX;
	b.length_of_longest_seq  = 0;
	std::cerr << "   ";
	while( mem < params->memory_limit && db_reader->has_next() ){
//		std::cerr << "\nmem: " << mem << "\n";
		std::cerr << ".";
		size_t c=0;
		while( c < params->sequence_block_size && db_reader->has_next() ){
			// s loaded with matrix
			Sequence *s = db_reader->get_next(false, false);
			++c;
			if( s->length()/shortest >= params->coverage_seq_len ){
				if (! params->profile_query) _6mer->create_kmer_list_fast_spaced( s );
				else _6mer->create_kmer_list_spaced(s);
				matches_6mer += table->match( _6mer );

				Base_table::MatchIterator it  = table->begin( (float)s->length(),
						params->kmer_score_threshold,
						params->coverage_seq_len );
				Base_table::MatchIterator end = table->end();
				if( it!=end ){
					_4mer_table->reset();
					if (! params->profile_query) _4mer->create_kmer_list_fast_spaced( s );
					else _4mer->create_kmer_list_spaced(s);

					size_t count_homologes=0;
					while( it!=table->end() && count_homologes<params->max_number_of_4mer_table_representants ){
						Sequence *r = block_man->get_representative( it.get_sequence_index() );
						_4mer_table->add_representative_spaced( r );
						++it;
						++count_homologes;
					}

					if( count_homologes>counter__max_kmeraln_count ){
						counter__max_kmeraln_count = count_homologes;
					}

					if( it!=table->end() && count_homologes==params->max_number_of_4mer_table_representants )
						throw MyException("Sequence limit of 4-mer table reached! (%i) ", params->max_number_of_4mer_table_representants );


					matches_4mer += _4mer_table->match_full( _4mer );
					Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
					Recycle_table::MatchIterator alns_end = _4mer_table->end();
					size_t idx=0;
					float score=0.0f;
					while( alns_it!=alns_end ){
						Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );
						_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r);
						float current_score = _4meraln->get_rel_score();
						float current_evalue = _4meraln->get_evalue();
						//						if(  current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_short()>=params->aln_coverage_short && _4meraln->get_aln_cov_long()>=params->aln_coverage_long ){
						if(  current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_long()>=params->aln_coverage_long ){
							score = current_score;
							idx = alns_it.get_sequence_index();
							break;
						}
						++alns_it;
					}
					if( idx>0 ){
						clusters->add_member( idx, s->index(), score );
						delete s;
						continue;
					}
				}
			}
			//temporary values for all sequences in this block
			//->are corrected in the final selfcompare method to take values of representantive sequences only
			b.length_of_shortest_seq = std::min( b.length_of_shortest_seq, s->length() );
			b.length_of_longest_seq  = std::max( b.length_of_longest_seq, s->length() );
			chars += s->length();
			list_of_new_candidates.push_back( s );
			seq_mem += s->get_memory_usage() + 2 * sizeof(void*);
		}
		mem = offset_memory_usage + block_man->get_memory_usage() + chars*sizeof( Base_table::_simple_rkmer ) + seq_mem;
//		std::cout << "Memory usage after while in _compare_db_reader: " << mem << ", seq_mem: " << seq_mem << "\n";
	}
//	std::cout << "Memory full: mem " << mem << ", memory limit: " << params->memory_limit << "\n";
	std::cerr << std::endl;
	b.characters = chars;
	b.end        = db_reader->get_index_of_next_sequence();
}

void K_Clust::_compare_list(){
	std::cout << "Function: _compare_list\n";
	std::list<Sequence*>::iterator cand_it  = list_of_new_candidates.begin();
	std::list<Sequence*>::iterator cand_end = list_of_new_candidates.end();
	const double shortest = block_man->get_current_table_block().length_of_shortest_seq;
	size_t c=0;
	std::cerr << "   ";
	while( cand_it!=cand_end ){
		++c;
		if( c%params->sequence_block_size==0 ) std::cerr << ".";
		Sequence *s = (*cand_it);

		if( s->length()/shortest < params->coverage_seq_len ){
			std::cerr << "   Remaining sequences are not in coverage range!" << std::endl;
			break;
		}

		if (! params->profile_query) _6mer->create_kmer_list_fast_spaced( s );
		else _6mer->create_kmer_list_spaced(s);
		matches_6mer += table->match( _6mer );

		Base_table::MatchIterator it  = table->begin( (float)s->length(), params->kmer_score_threshold, params->coverage_seq_len );
		Base_table::MatchIterator end = table->end();
		if( it!=end ){
			_4mer_table->reset();
			if (! params->profile_query) _4mer->create_kmer_list_fast_spaced( s );
			else _4mer->create_kmer_list_spaced(s);

			size_t count_homologes=0;
			while( it!=table->end() && count_homologes<params->max_number_of_4mer_table_representants ){
				Sequence *r = block_man->get_representative( it.get_sequence_index() );
				_4mer_table->add_representative_spaced( r );
				++it;
				++count_homologes;
			}

			if( count_homologes>counter__max_kmeraln_count ){
				counter__max_kmeraln_count = count_homologes;
			}

			if( it!=table->end() && count_homologes==params->max_number_of_4mer_table_representants )
				throw MyException("Sequence limit of 4-mer table reached! (%i) ", params->max_number_of_4mer_table_representants);

			matches_4mer += _4mer_table->match_full( _4mer );
			Recycle_table::MatchIterator alns_it  = _4mer_table->begin();
			Recycle_table::MatchIterator alns_end = _4mer_table->end();
			size_t idx=0;
			float score=0.0f;
			while( alns_it!=alns_end ){
				Sequence *r = block_man->get_representative( alns_it.get_sequence_index() );

				_4meraln->align( Kmeraln::FAST_ADDR, alns_it.get_match_list(), s, r);
				float current_score = _4meraln->get_rel_score();
				float current_evalue = _4meraln->get_evalue();

				//				if( current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_short()>=params->aln_coverage_short && _4meraln->get_aln_cov_long()>=params->aln_coverage_long){
				if( current_score>params->clustering_threshold && current_score>score && current_evalue<=params->clustering_evalue_threshold && _4meraln->get_aln_cov_long()>=params->aln_coverage_long){
					score = current_score;
					idx = alns_it.get_sequence_index();
					break;
				}
				++alns_it;
			}
			if( idx>0 ){
				clusters->add_member( idx, s->index(), score );
				delete s;
				std::list<Sequence*>::iterator tmp = cand_it;
				++cand_it;
				list_of_new_candidates.erase(tmp);
				continue;
			}
		}
		++cand_it;
	}
	std::cerr << std::endl;
}

void K_Clust::free_tables_seedlists_kmeraln(){

	//std::cerr << "delete table" << std::endl;
	delete table;
	//std::cerr << "delete _6mer" << std::endl;
	delete _6mer;
	//std::cerr << "delete _4mer" << std::endl;
	delete _4mer;
	//std::cerr << "delete _4meraln" << std::endl;
	//delete _4meraln;
	//std::cerr << "delete _4mer_table" << std::endl;
	delete _4mer_table;

	//std::cerr << "set all = 0" << std::endl;
	table       = 0;
	_6mer       = 0;
	_4mer       = 0;
	_4meraln    = 0;
	_4mer_table = 0;

}

void K_Clust::cleanup_tmp_files() throw (std::exception){
	std::cerr << "======================================================" << std::endl;
	if( params->remove_table_files && (block_man->get_number_of_blocks()!=1) ){
		BlockManager::Iterator it = block_man->begin();
		while( it!=block_man->end() ){
			std::cerr << "   deleting " << (*it).get_filename() << std::endl;
			if( remove((*it).get_filename().c_str()) )
				throw MyException( 1, "Cannot delete '%s'\n", (*it).get_filename().c_str() );
			++it;
		}
	}

	if( params->remove_dbsorted ){
		std::cerr << "   deleting " << params->get_dbsortedfile() << std::endl;
		if( remove(params->get_dbsortedfile().c_str()) ){
			throw MyException( 1, "Cannot delete '%s'\n", params->get_dbsortedfile().c_str() );
		}
	}
	std::cerr << std::endl;
}

void K_Clust::write_representatives_db() throw (std::exception){
	clusters->write_representants_db( db_reader,
			params->get_representants_dbfile().c_str(),
			params->representatives_header_merging_method);
}

void K_Clust::write_dmp_files() throw (std::exception){
	clusters->write_dmp_files( db_reader,
			params->get_header_dmp_file().c_str(),
			params->get_nodes_dmp_file().c_str(),
			params->profile_query);
}

void K_Clust::write_refined_nodes_dmp() throw (std::exception){
	if (!params->profile_query)
		clusters->write_refined_nodes_dmp_file( params->get_refined_nodes_dmp_file().c_str());
	else
		clusters->write_dmp_files_profile(db_reader, params->get_refined_nodes_dmp_file().c_str());
}

void K_Clust::write_debug_files(const char *fn) throw (std::exception){
	//debug
	//std::string fn = params->get_working_dir() + "repnumbers.dat";
	//clusters->write_representative_numbers( fn.c_str() );
	clusters->write_debug_file( db_reader, fn);
}

void K_Clust::insert_into_queue(std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare>* q, time_benchmark* t){
	if (t->time < 500) return;
	q->push(t);
	if (q->size() > benchmark_queue_size){
		time_benchmark* tt = q->top();
		q->pop();
		delete tt;
	}
}

void K_Clust::print_memory_usage()
{
	   using std::ios_base;
	   using std::ifstream;
	   using std::string;

	   double vm_usage     = 0.0;
	   double resident_set = 0.0;

	   // 'file' stat seems to give the most reliable results
	   ifstream stat_stream("/proc/self/stat",ios_base::in);

	   // dummy vars for leading entries in stat that we don't care about
	   //
	   string pid, comm, state, ppid, pgrp, session, tty_nr;
	   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	   string utime, stime, cutime, cstime, priority, nice;
	   string O, itrealvalue, starttime;

	   // the two fields we want
	   //
	   unsigned long vsize;
	   long rss;

	   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	               >> utime >> stime >> cutime >> cstime >> priority >> nice
	               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	   stat_stream.close();

	   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	   vm_usage     = vsize / 1024.0;
	   resident_set = rss * page_size_kb;

	   std::cout << "--------------------------------------\n";
	   std::cout << "VM      : " << vm_usage << "\n";
	   std::cout << "RESIDENT: " << resident_set << "\n";
	   std::cout << "KCLUST  : " << get_memory_usage() << "\n";
	   std::cout << "--------------------------------------\n";

}


