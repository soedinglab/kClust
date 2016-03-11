/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include <vector>
#include <iostream>
#include <string>

#include "mysqldb.h"
#include "params.h"

using std::endl;

std::ostream& usage( std::ostream &out );
void set_params( const int, const char**, Params& ) throw (std::exception);

int main(int argc, const char **argv){
	if( argc < 5 ){
		usage(std::cerr);
		return 1;
	}

	try{
		Params params;
		set_params(argc, argv, params);

		Matrix m(params.matrix_type);

		std::cerr << "Init..." << std::endl;

		std::cerr << "   MySQL ..." << std::endl;
		Mysqldb mdb( &params, &m);
		std::cerr << "ok" << std::endl;

		std::cerr << "   Database reader...";
		Fasta_db_reader dbr( params.get_dbsortedfile().c_str(), &m);
		std::cerr << "ok" << std::endl;

		std::cerr << "   Computing alignments..." << std::endl;
		mdb.createdb( &dbr, params.consensus_header_merging_method );
		std::cerr << "   done." << std::endl;

		std::cerr << "   Writing consensus db...";
     	mdb.write_consensus_db_file();
      std::cerr << "   ok" << std::endl;
		
		std::cerr << "done." << std::endl;

	}catch(std::exception &e){
		std::cerr << std::endl << ">" << std::endl;
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << ">" << std::endl;
		return 2;
	}

  	return 0;
}

void set_params(const int argc, const char **argv, Params &params) throw (std::exception){
	for(int i=1; i<argc; ++i){
		std::string arg(argv[i]);
		if( arg=="-c" ){
			if( ++i<argc ) params.alignment_cmd = std::string(argv[i]);
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-d" ){
			if( ++i<argc ){ 
				params.set_working_dir( argv[i] );
				if( !Params::exists( params.get_dbsortedfile().c_str() ) ) 
					throw MyException("Cannot find database file: '%s' in directory: '%s'", params.get_dbsortedfile().c_str(), params.get_working_dir().c_str());
				if( Params::exists( params.get_refined_nodes_dmp_file().c_str() ) )
					std::cerr << "File with refined clustering found: " << params.get_refined_nodes_dmp_file() << std::endl;
				else if( !Params::exists( params.get_nodes_dmp_file().c_str() ) ) 
					throw MyException("Cannot find cluster file '%s' or '%s'", params.get_refined_nodes_dmp_file().c_str(), params.get_nodes_dmp_file().c_str());
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="--mysql-host" ){
			if( ++i<argc ) params.mysql_host = std::string(argv[i]);
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--mysql-db" ){
			if( ++i<argc ) params.mysql_db = std::string(argv[i]);
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--mysql-user" ){
			if( ++i<argc ) params.mysql_user = std::string(argv[i]);
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--mysql-pass" ){
			if( ++i<argc ) params.mysql_passwd = std::string(argv[i]);
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="-r" ){
			params.continue_aln_comp = true;
			continue;
		}

		if( arg=="--merge-ncbi-headers" ){
			params.consensus_header_merging_method = Clusters::NCBI_HEADER_MERGING;
			continue;
		}
		
		if( arg=="--no-pseudo-headers" ){
			params.write_pseudo_headers = false;
			continue;
		}

		if( arg=="--a3m" ){
			params.alignment_output_format = Alignment::a3m;
			continue;
		}

		if( arg=="-h" || arg=="-help" || arg=="--h" || arg=="--help" ){ 
			usage(std::cerr);
			exit(2);
		}
		throw MyException("Unknown command-line parameter found: '%s'", argv[i]);
	}
}

std::ostream& usage( std::ostream &out ){

	char buffer[512];
	sprintf(buffer, "Usage: ./kClust_mkAlnDB -c '[command]' -d [directory] [options]\n");
	out << buffer;
	out << std::endl;

	out << "Version 1.0" << endl << endl;

	out << "kClust_mkAlnDB creates an alignment database." << endl;
	out << "Written by Christian Mayer (christian.eberhard.mayer@googlemail.com)" << endl << endl;


	out << "Required arguments:" << endl;
	sprintf(buffer, " %-20s %-20s: %s\n", "-c", "'[command]'", "Command for multiple sequence alignment program (input/output must be in FastA format)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-d", "[directory]", "Directory with the result files of kClust.");
	out << buffer;
	out << std::endl;
	out << "Optional arguments:" << endl;
	sprintf(buffer," %-20s %-20s: %s\n", "-r", "", "Recover program abortion.");
	out << buffer;
	sprintf(buffer," %-20s %-20s: %s\n", "--mysql-host", "[host]", "Hostname where the MySQL server resides (default=localhost).");
	out << buffer;
	sprintf(buffer," %-20s %-20s: %s\n", "--mysql-db", "[name]", "Name of MySQL database (default=nr30).");
	out << buffer;
	sprintf(buffer," %-20s %-20s: %s\n", "--mysql-user", "[user]", "Username for connection to MySQL server (default=nr30).");
	out << buffer;
	sprintf(buffer," %-20s %-20s: %s\n", "--mysql-pass", "[password]", "Password for connection to MySQL server (default=nr30).");
	out << buffer;
	sprintf(buffer," %-20s %-20s: %s\n", "--a3m", "", "Save alignments in a3m-format (default=fasta).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--merge-ncbi-headers", "", "Compress NCBI headers for consensus sequences (default=concatenation).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--no-pseudo-headers", "", "Write original header in cluster alignments (default=pseudo-headers).");
	out << buffer;
	out << std::endl;
	out << " Example: ./kClust_mkAlnDB -c '/usr/bin/kalign -q -i $infile -o $outfile' -d /tmp/clustering/" << endl << endl;
	return out;
}
