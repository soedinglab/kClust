#ifndef CM_CLUSTERS_H
#define CM_CLUSTERS_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Maintains the clusters of representatives

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
//#include <algorithm>

#include "simple_list.h"
#include "simple_hash.h"
#include "fasta_db_reader.h"
#include "converter.h"
#include "profile.h"

class Clusters{
	
	friend class Mysqldb;
	friend class MkAln;

	public:

		typedef std::vector<std::string> string_vector_type;
		class Accessions{
			public:
				Accessions():name(""),prio(0){}
				static bool greater(const Accessions &a, const Accessions &b){ return a.prio>b.prio; }
				std::string name;
				string_vector_type gis;
				string_vector_type pdbs;
				string_vector_type sps;
				int prio;
		};	

		class Member{
			public:
				Member( const size_t idx=0, const float score=0.0f ):idx(idx), score(score){}
				size_t idx;
				float score;
		};

		enum header_merging_method_t{
			SIMPLE_CONCATENATION=0,
			NCBI_HEADER_MERGING=1,
            UNIPROT_HEADER_MERGING=2,
            ORIGINAL_REPRESENTATIVE_HEADER=3
		};

		class prio_string{
			public:
				prio_string(std::string s="", int p=0):str(s),prio(p){}
				prio_string(const prio_string &other):str(other.str),prio(other.prio){}
				bool operator<(const prio_string &other)const{ return prio<other.prio;}
				static bool greater(const prio_string &a, const prio_string &b){ return a.prio>b.prio; }
				bool operator==(const prio_string &other)const{ 
					return prio==other.prio && str==other.str; 
				}
				bool operator!=(const prio_string &other)const{ 
					return prio!=other.prio || str!=other.str; 
				}
				prio_string& operator=(const prio_string &other){ 
					str = other.str;
					prio = other.prio;
					return *this;
				}

				std::string str;
				int prio;
		};


		Clusters(const size_t dbsize, Matrix *m);
		Clusters(const char *clusters_dmp_file, Matrix *m);
		virtual ~Clusters();

		bool is_representative_in( const size_t seqindex, const size_t blockindex );

		bool is_representative( const size_t idx );

		// idx_r: cluster index, idx_m: sequence index, score: score of the sequence idx_m with the representative of the cluster
		void add_member( const size_t idx_r, const size_t idx_m, const float score );
		
		// index: index of the sequence, blockindex: index of the block that contains the sequence
		void create_cluster_for( const size_t index, const size_t blockindex );

		void evaluate_members();

		void remove_members();

		void print(Fasta_db_reader *dbr);

		void write_representants_db( 	Fasta_db_reader *dbr, 
						const char *outfile, 
						const header_merging_method_t merging_method) throw (std::exception);

		void write_dmp_files(	Fasta_db_reader *dbr, 
					const char *names, 
					const char *clusters,
					const bool profile_query) throw (std::exception);

		void write_dmp_files( Fasta_db_reader *dbr, const char *clusters, const char *names = NULL) throw (std::exception);

		void write_refined_nodes_dmp_file(const char *clusters) throw (std::exception);
		
		void write_dmp_files_profile(Fasta_db_reader *dbr, const char *clusters, const char *headers = NULL) throw (std::exception);


		//numbers of representatives - swapping debug
		void write_representative_numbers(const char *filename) throw (std::exception);

		//prints representative headers and member headers
		void write_debug_file(Fasta_db_reader *dbr, const char *filename) throw (std::exception);

		size_t get_number_of_clusters();

		size_t get_member_count( const size_t idx );

		void make_aln(const char *infile, const char *outfile);

		const size_t get_memory_usage() const;

		std::ostream& print_representants_with_members( std::ostream &out, Fasta_db_reader *dbr );
		
//		void analyze_keywords(Fasta_db_reader *dbr) throw(std::exception);
		
	private:
		size_t dbsize;
		Matrix *matrix;
		Simple_list<Member> **representatives;
		size_t *belongs_to_block;
		size_t _cluster_count;
	
		void eval_ncbi_header(	const char *header,
					std::vector<Accessions> &org_vec,
					Simple_hash<size_t> &org_check,
					Simple_hash<prio_string> &check
				) throw (std::exception);

		void eval_uniprot_header(	const char *header,
							std::vector<Accessions> &org_vec,
							Simple_hash<size_t> &org_check,
							Simple_hash<prio_string> &check
						) throw (std::exception);

		void _init_mysql_tables() throw (std::exception);

		std::ostream& create_merged_header_simple( std::ostream &hstream, string_vector_type &headers, Sequence* r, const size_t cluster_size);

		std::ostream& create_merged_header_advanced( std::ostream &hstream, 
					  		     Sequence *r, 
							     const size_t cluster_size,
							     std::vector<Accessions> &org_vec,
							     Simple_hash<size_t> &org_check,
							     Simple_hash<prio_string> &check );

		void read_clusters_dmp_file(const char *filename ) throw (std::exception);
		void write_dmp_files_np(Fasta_db_reader *dbr, const char *clusters, const char *names) throw (std::exception);
		
};

#endif
