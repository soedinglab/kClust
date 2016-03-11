/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include "params.h"
#include "my_exception.h"

Params::Params(){
	matrix_type                            = Matrix::static_blosum62;
	db                                     = "";
	input_dir                              = "";
	_4mer_m_file                            = "";
	_2mer_m_file                            = "";
	
	hash_bits                              = 26;
	memory_limit                           = 1024u*1024u*3500u; // 1408987326;
	max_number_of_4mer_table_representants = 100000;
//	max_number_of_4mer_table_representants = 10000;
	kmer_score_threshold                   = 0.55f; //98%
	clustering_threshold                   = 1.12f;
	clustering_evalue_threshold            = 1.00e-3f;
	kdp_kmer_similarity_threshold          = 2.9f;
//	kdp_kmer_similarity_threshold          = 2.8f;
	kdp_k                                  = 4;
	kdp_p_m                                = 1.00e-4f;
	kdp_G                                  = 12.0;
	kdp_E                                  = 2.0;
	kdp_F                                  = 0.27;
	kdp_delta                              = 50;
	kdp_Lamda                              = 0.3;
	kdp_K                                  = 0.09;
	kdp_H                                  = 0.34;
	filter_k                               = 6;
	filter_kmer_similarity_threshold       = 4.3f;
//	filter_kmer_similarity_threshold       = 4.2f;
	check_coverage_seq_len                 = true;
	coverage_seq_len                       = 0.8f;
	aln_coverage_long                      = 0.8f;
	remove_dbsorted                        = false;
	remove_table_files                     = false;
	representatives_header_merging_method  = Clusters::ORIGINAL_REPRESENTATIVE_HEADER;
	consensus_header_merging_method        = Clusters::SIMPLE_CONCATENATION;
	refinement                             = false;
	profile_query                          = false;
	score_correction                       = false;
	kmer_matrices                          = false;
	sequence_block_size                    = 10000;
	db_length                              = 0;
	sum_log_db_length                      = 0;
	
	write_add_header					   = false;
	write_pseudo_headers				   = true;
	advanced_clustering                    = true;
	
	mysql_host                             = "localhost";
	mysql_user                             = "nr30";
	mysql_passwd                           = "nr30";
	mysql_db                               = "nr30";
	mysql_seq_table_name                   = "sequences";
	mysql_clu_table_name                   = "clusters";

	continue_aln_comp                      = false;

	alignment_output_format                = Alignment::fasta;

	consensus_method					   = 1;

	working_dir                            = "";
	tmp_dir                                = "";
	aln_db_dir 							   = "";
	aln_dir_old                            = "";
	
	cons_db_file                           = "";
	parallel_aln_comp                      = false;
	
	benchmark                              = -1;
	count_idents                           = false;

	kClust_cmd                             = "/cluster/toolkit/production/bioprogs/kClust/bin/kClust";

	write_time_bench                       = false;
}

void Params::set_working_dir(const std::string &dir) throw (std::exception){
	working_dir = dir;
	if( working_dir.at(working_dir.length()-1 )!='/' ) working_dir += "/";
	if( !Params::is_directory( working_dir.c_str() ) ){
		std::string cmd = "mkdir " + working_dir;
		int i = system( cmd.c_str() );
		if (i != 0) throw MyException("The working directory '%s' could not be created!", working_dir.c_str());
		//throw MyException("The denoted working directory '%s' does not exist or is not a directory!", working_dir.c_str());
	}
	
	if (profile_query) dbsorted = working_dir + "db_cons_sorted.fas";
	else dbsorted = working_dir + "db_sorted.fas";
	
	if (profile_query) db = working_dir + "db_tmp.fas";
	header_dmp_file         = working_dir + "headers.dmp";
	nodes_dmp_file          = working_dir + "clusters.dmp";
	refined_nodes_dmp_file  = working_dir + "refined_clusters.dmp";
	rep_db_file             = working_dir + "representatives.fas";
	if (cons_db_file == "") {
		cons_db_file        = working_dir + "consensus.fas";
	}
	set_tmp_dir(working_dir + "tmp" + "/");
}

void Params::set_tmp_dir(const std::string &tmpd) throw (std::exception){
	if (tmp_dir != ""){
		std::string cmd = "rm -r " + tmp_dir;
		system( cmd.c_str() );
	}
	tmp_dir = tmpd;
	if (!Params::is_directory(tmp_dir.c_str())){
		std::string cmd = "mkdir " + tmp_dir;
		system( cmd.c_str() );
	}
}

void Params::set_input_fileordir(const std::string &fileordir) throw (std::exception){
	if (profile_query){
		Params::set_input_dir(fileordir);
	}
	else{
		Params::set_dbfile(fileordir);
	}
}

void Params::set_aln_db_dir(const std::string &dir) throw (std::exception){
	aln_db_dir = dir;
	
	if( aln_db_dir.at(aln_db_dir.length()-1 )!='/' ) aln_db_dir += "/";
	
	if( Params::is_directory( aln_db_dir.c_str() ) ) { 
		//throw MyException("'%s' does already exist!", aln_db_dir.c_str());
		std::cerr << "WARNING! " << aln_db_dir << " already exists!" << std::endl;
		return;
	} else {
		std::string cmd = "mkdir " + aln_db_dir;
		int ret = system( cmd.c_str() );
		if (ret != 0)
			throw MyException("Cannot create '%s'!", aln_db_dir.c_str());
	}
	// Create subdirectories
	char first[] = {'B', 'C', 'D', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z'};
	char second[] = {'A', 'E', 'I', 'O', 'U'};
	
	for (unsigned int i = 0; i < 20; i++) {
		for (unsigned int j = 0; j < 5; j++)  {
			std::string cmd = "mkdir " + aln_db_dir + first[i] + second[j];
			system( cmd.c_str() );
		}
	} 
	
		std::cout << "setting aln_db_dir to " << aln_db_dir << "\n";
	
}	

std::string Params::get_working_dir(){ return working_dir; }

std::string Params::get_aln_db_dir(){ 
	return aln_db_dir; 
	}

void Params::set_consensus_dbfile(const std::string &dbfile) { cons_db_file = dbfile; }

void Params::set_dbfile(const std::string &dbfile) throw (std::exception){
	db = dbfile;
	if( !Params::exists( db.c_str() ) ) 
		throw MyException("The database file '%s' does not exist!", db.c_str());
}

void Params::set_input_dir(const std::string &input_dir_) throw (std::exception){
	input_dir = input_dir_;
	if( input_dir.at(input_dir.length()-1 )!='/' ) input_dir += "/";
	
	clusters_old = input_dir + "clusters.dmp";
	headers_old = input_dir + "headers.dmp";
	aln_dir_old = input_dir + "alignments/";
	dbsorted_old = input_dir + "db_sorted.fas";
		
	if( !Params::exists( clusters_old.c_str() ) ) 
		throw MyException("The clusters file '%s' in the input directory does not exist!", clusters_old.c_str());
	
	if( !Params::exists( aln_dir_old.c_str() ) ) 
		throw MyException("No directory with alignments '%s' exists in the input directory!", aln_dir_old.c_str());
	
	if( !Params::exists( headers_old.c_str() ) ) 
		throw MyException("The headers file '%s' in the input directory does not exist!", headers_old.c_str());
	
//	if( !Params::exists( dbsorted_old.c_str() ) ) 
//		throw MyException("The sorted db file '%s' in the input directory does not exist!", headers_old.c_str());
	
}

std::string Params::get_dbfile(){ return db; }

std::string Params::get_input_dir(){ return input_dir; }

std::string Params::get_tmp_dir(){ return tmp_dir; }

std::string Params::get_dbsortedfile(){ return dbsorted; }

std::string Params::get_clustersfile_old(){ return clusters_old; }

std::string Params::get_dbsortedfile_old(){ return dbsorted_old; }

std::string Params::get_headersfile_old(){ return headers_old; }

std::string Params::get_aln_dir_old(){ return aln_dir_old; }

std::string Params::get_header_dmp_file(){ return header_dmp_file; }

std::string Params::get_nodes_dmp_file(){ return nodes_dmp_file; }

std::string Params::get_refined_nodes_dmp_file(){ return refined_nodes_dmp_file; }

std::string Params::get_representants_dbfile(){ return rep_db_file; }

std::string Params::get_consensus_dbfile(){ return cons_db_file; }

bool Params::exists(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		ret=true;
	return ret;
}

bool Params::has_zero_size(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		if(fstats.st_size==0) ret=true;
	return ret;
}

bool Params::is_directory(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		if( S_ISDIR(fstats.st_mode)) ret=true;
	
	return ret;
}

bool Params::is_link(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		if( S_ISLNK(fstats.st_mode)) ret=true;
	
	return ret;
}

bool Params::is_regular_file(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		if( S_ISREG(fstats.st_mode)) ret=true;
	
	return ret;
}

void Params::check() throw (std::exception){
	if( working_dir == "" ) throw MyException("Directory for results and tmp-files not set.");
	if( db == "" && profile_query == false)	throw MyException("No database specified.");
	if ( input_dir == "" && profile_query == true) throw MyException("No input directory specified.");
//	if ( profile_query == true && !Params::exists(clusters_old.c_str())) throw MyException("Input directory does not contain a clusters.dmp file.");
	if( kdp_k<1 || kdp_k>6)  throw MyException("k-mer length for kDP out of range [2..6]!");
	if( filter_k<1 || filter_k>6)  throw MyException("k-mer length for similarity score filter out of range [2..6]!");
	if (_4mer_m_file == "" && _2mer_m_file != "" || _4mer_m_file != "" && _2mer_m_file == "") throw MyException("Only one k-mer substitution matrix file specified!");
//	if (score_correction && profile_query) throw MyException("Sequence background frequency score correction for the k-mer scores not implemented for the profile queries.");
}

