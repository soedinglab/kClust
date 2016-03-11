#ifndef CM_MYSQLDB_H
#define CM_MYSQLDB_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION
//Implements the setup of an alignment and consensus sequence database.
//Creates two relational MySQL tables.

#include <vector>
#include <string>
#include "clusters.h"
#include "params.h"

//forward declarations - MySQL header files cannot be included here
struct mysql_cluster;
struct mysql_sequence;

class Mysqldb{

	public:
		Mysqldb(Params *params, Matrix *m) throw (std::exception);
		~Mysqldb();
		void createdb(	Fasta_db_reader *dbr, 
							const Clusters::header_merging_method_t merging_method) throw (std::exception);
		void write_consensus_db_file() throw (std::exception);
	private:
		Clusters clu;
		Params *params;
		std::string cmd_first_part;
		std::string cmd_second_part;
		std::string cmd_third_part;
		size_t start_index;
		size_t cl_idx;
		void init_mysql_tables() throw (std::exception);
		void recover() throw (std::exception);
		void insert_seq_into_db( mysql_sequence &mseq_struct ) throw (std::exception);
		void insert_cluster_into_db( mysql_cluster &cl ) throw (std::exception);
};

#endif
