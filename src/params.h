#ifndef CM_PARAMS_H
#define CM_PARAMS_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Parameter wrapper

#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "matrix.h"
#include "clusters.h"
#include "alignment.h"

class Params{

	public:
		Params();
		void check() throw (std::exception);
		void set_working_dir(const std::string &dir) throw (std::exception);
		std::string get_working_dir();
		void set_input_fileordir(const std::string &fileordir) throw (std::exception);
		void set_dbfile(const std::string &dbfile) throw (std::exception);
		void set_input_dir(const std::string &input_dir_) throw (std::exception);
		void set_tmp_dir(const std::string &tmp_dir_) throw (std::exception);
		std::string get_dbfile();
		std::string get_input_dir();
		
		void set_aln_db_dir(const std::string &dir) throw (std::exception);
		std::string get_aln_db_dir();
		
		void set_consensus_dbfile(const std::string &dbfile);
		
		std::string get_dbsortedfile();
		std::string get_clustersfile_old();
		std::string get_headersfile_old();
		std::string get_aln_dir_old();
		std::string get_dbsortedfile_old();
		std::string get_representants_dbfile();
		std::string get_consensus_dbfile();
		std::string get_header_dmp_file();
		std::string get_nodes_dmp_file();
		std::string get_refined_nodes_dmp_file();
		std::string get_tmp_dir();

		static bool has_zero_size(const char *);
		static bool is_directory(const char *);
		static bool is_regular_file(const char *);
		static bool is_link(const char *);
		static bool exists(const char *);

		//type of matrix to be used e.g. blosum62
		Matrix::mtype matrix_type;

		//bit-length of hash-keys for clustering at 100%
		size_t hash_bits;

		//maximum number of bytes that can be used by the program
		size_t memory_limit;

		//length of k-mers for kDP
		size_t kdp_k;
		
		//kdp delta window 
		int kdp_delta;

		//probability for chance match
		float kdp_p_m;
	
		//gap open penalty
		float kdp_G;

		//gap extension penalty
		float kdp_E;
	
		//intra-diagonal gap penalty
		float kdp_F;
		
		// Statistical Parameters for BLOSUM62 matrix
		float kdp_Lamda;
		float kdp_K;
		float kdp_H;

		//length of k-mers for k-mer similarity scoring filter
		size_t filter_k;

		//threshold 6mer-score/length-of-query -> triggers kmer alignments
		float kmer_score_threshold;

		//threshold for seed generation of 6mers
		float filter_kmer_similarity_threshold;

		//threshold for generation of kmers for kDP alignments
		float kdp_kmer_similarity_threshold;

		//minimum sequence identity in bits per column of the shorter sequence 
		float clustering_threshold;
		
		//maximum E-value for clustering 
		float clustering_evalue_threshold;

		//triggers optimal assignment of members to representatives
		bool refinement;
				
		// the query is a profile and the corresponding consensus sequence (instead of a simple sequence)
		bool profile_query;
		
		bool score_correction;
		
		// use 4-mer and 2-mer matrices
		bool kmer_matrices;
				
		//4mer alignments limit
		//maximum number of representants that can be hold by the 4mer table
		size_t max_number_of_4mer_table_representants;

		//coverage
		bool check_coverage_seq_len;
		float coverage_seq_len;

		//the fraction of residues of the longer sequence that have to be aligned
		float aln_coverage_long;
		
		//the fraction of residues of the shorter sequence that have to be aligned
//		float aln_coverage_short;

		//output alignment format in MySQL db 
		Alignment::format alignment_output_format;

		//consensus calculation method
		// 0: no consensus sequence
		// 1: conventional method (orig. code)
		// 2: max. profile score with consensus
		// 3: min. rel. entropy between p(a) and P(a|x)
		// 4: min. quadratic diff. between profile score and matrix score
		// 5: min. expected score difference between homologs and random seqs
		int consensus_method;
		
		//sequence reading block size
		size_t sequence_block_size;

		//flag if files should be deleted
		bool remove_table_files;
		bool remove_dbsorted;
		
		//method for header merging
		Clusters::header_merging_method_t representatives_header_merging_method;
		
		//method for header merging
		Clusters::header_merging_method_t consensus_header_merging_method;

		//alignment program
		std::string alignment_cmd;
		
		//kClust location
		std::string kClust_cmd;

		//kClust_mkAlb location
		std::string kClust_mkAln_cmd;

		// write_pseudo_headers
		bool write_pseudo_headers;

		// write additional header in alignment
		bool write_add_header;

		// perform advanced clustering (recluster cluster with more than 8000 sequences with 90% seq-ID)
		bool advanced_clustering;

		std::string mysql_host;
		std::string mysql_user;
		std::string mysql_db;
		std::string mysql_passwd;
		std::string mysql_seq_table_name;
		std::string mysql_clu_table_name;

		bool continue_aln_comp;
		
		// database length (L_DB)
		float db_length;
		
		// sum of log (template-length) - needed for DB-length correction   (1/H SUM log(L_T))
		float sum_log_db_length;
		
		//directory with alignments and sorted database of the previous kClust run for profile queries
		std::string input_dir;
		
		//============== only for benchmarks ========================
		// benchmark type
		int benchmark;
		
		// count identities instead of similar k-mers
		bool count_idents;

		std::string _4mer_m_file;

		std::string _2mer_m_file;

		// parallelized alignment computation on the cluster
		bool parallel_aln_comp;

		bool write_time_bench;

	private:
		//filename of working-directory
		std::string working_dir;

		//filename of directory with alignment files
		std::string aln_db_dir;

		//filename of fasta database
		std::string db;

		//directory with temporary data (swapped tables)
		std::string tmp_dir;

		//filename of fasta database
		std::string dbsorted;
				
		//for profile query, path of the headers.dmp file of the previous kClust run
		std::string headers_old;
		
		//for profile query, path of clusters file of the previous kClust run
		std::string clusters_old;
		
		//for profile query, path of the dir with alignments of the previous kClust run
		std::string aln_dir_old;
		
		//for profile query, path to the sorted db of the previous kClust run
		std::string dbsorted_old;
		
		//filename of representatives database
		std::string rep_db_file;
		
		//filename of consensus db
		std::string cons_db_file;

		//filename of headers dump
		std::string header_dmp_file;

		//filename of nodes dump
		std::string nodes_dmp_file;

		//filename of refined nodes dump
		std::string refined_nodes_dmp_file;
		
};

#endif
