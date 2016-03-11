#ifndef CM_KCLUST_H
#define CM_KCLUST_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//wrapper class combining the k-mer similarity scoring and
//k-mer dynamic programming method
//into the kClust database clustering method.
//swapping is employed if required.

#include<string>
using std::string;

#include <list>
#include <vector>
#include <iostream>
#include <iomanip>
#include <queue>

#include "matrix.h"
#include "fasta_db_reader.h"
#include "base_table.h"
#include "block_manager.h"
#include "clusters.h"
#include "kmeraln.h"
#include "recycle_table.h"
#include "params.h"
#include "sequence.h"
#include "converter.h"
#include "chronometer.h"

//#define __CM_DEBUG_PRINT_LETTERS 1

class K_Clust{
	public:

		class _edge{
			public:
				_edge(int source=0, int target=0):source(source), target(target){};
				int source;
				int target;
		};
		
		struct time_benchmark{
			public:
			
			time_benchmark(int seq_id, int length, float _6mer_per_pos, float _4mer_per_pos, int aln_num, double time, double time_per_pos):
			seq_id(seq_id), length(length), _6mer_per_pos(_6mer_per_pos), _4mer_per_pos(_4mer_per_pos), aln_num(aln_num), time(time), time_per_pos(time_per_pos){}
			
			int seq_id; 
			int length;
			float _6mer_per_pos;
			float _4mer_per_pos;
			int aln_num;
			double time; 
			double time_per_pos;
		};
		
		struct time_benchmark_compare{
			bool operator() ( const time_benchmark* first, const time_benchmark* second ) const{
				return (first->time_per_pos > second->time_per_pos);
			}
		};

		K_Clust(Params*) throw (std::exception);
		~K_Clust();

		// create temporary consensus sequences database
		void prepare(const size_t width=80) throw (std::exception);
		
		void cluster() throw (std::exception);
		void refine() throw (std::exception);

		void write_dmp_files()  throw (std::exception);

		void write_refined_nodes_dmp()  throw (std::exception);
		void write_representatives_db()  throw (std::exception);
		void write_debug_files(const char *fn) throw (std::exception);

		void cleanup_tmp_files() throw (std::exception);

		void free_tables_seedlists_kmeraln();

		//std::ostream& write_graphml_file(std::ostream& out);
	
		size_t get_memory_usage();
		
		void pref_benchmark_scop_tpfp() throw (std::exception);
		void pref_benchmark_scop_scatterplot() throw (std::exception);
		void kDP_benchmark_scop_tpfp() throw (std::exception);
		void kDP_benchmark_scop_TMscore() throw (std::exception);

	private:
		Params *params;

		size_t block_index;
		size_t offset_memory_usage;
		size_t selfcompare_mem_limit;

		Fasta_db_reader *db_reader;
		BlockManager *block_man;
		Clusters *clusters;
		Base_table *table;
		Kmer *_6mer;
		Kmer *_4mer;
		Kmer *_3mer;
		Kmer_matrix* _4mer_matrix;
		Kmer_matrix* _2mer_matrix;
		Kmeraln *_4meraln;
		Recycle_table *_4mer_table;
		std::list<Sequence*> list_of_new_candidates;

		Matrix matrix;

		int node_count;
		
		void _init() throw (std::exception);
		// clustering of the block b
		void _selfcompare_db_reader(BlockManager::Block &b) throw (std::exception);
		// clustering of the block b
		void _selfcompare_list(BlockManager::Block &b) throw (std::exception);
		// Compares *s with all representative sequences
		// *s is added to the representatives if no representative is above kmer filter score
		// table_index: index of the block that contains the sequence
		void __selfcompare(	Sequence *s, 
									size_t &minlength, 
									size_t &maxlength, 
									const size_t table_index ) throw (std::exception);
		void _compare_db_reader(BlockManager::Block &b);
		void _compare_list();
		
		void insert_into_queue(std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare>* q, time_benchmark* t);
		
		void print_memory_usage();

		size_t counter__max_kmeraln_count;
		size_t counter__all_kmer_alns;
		size_t counter__kmeraln_trigger_seqs;
		
		double kmer6_av;
		double kmer4_av;
		size_t seqs_6;
		size_t seqs_4;
		long int matches_6mer;
		long int matches_4mer;
		
		int benchmark_queue_size;
		std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare>* q_pref;
		std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare>* q_kdp;
		std::priority_queue<time_benchmark*, std::vector<time_benchmark*>, time_benchmark_compare>* q_seq;

		
};
#endif
