/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 *                                                                         *
 *   rewritten by Michael Remmert (2008)                                   *
 *   remmert@lmb.uni-muenchen.de                                           *
 ***************************************************************************/
#include <vector>
#include <iostream>
#include <string>

#include "mkaln.h"
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

		MkAln mkaln( &params, &m);
		
		std::cerr << "Database reader...";
		Fasta_db_reader dbr( params.get_dbsortedfile().c_str(), &m, false, false);
		std::cerr << "ok" << std::endl;

		std::cerr << "Computing alignments..." << std::endl;
		if (params.parallel_aln_comp)
			mkaln.make_alignments_parallel( &dbr, params.consensus_header_merging_method );
		else
			mkaln.make_alignments( &dbr, params.consensus_header_merging_method );

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
		if (arg == "-cl"){
			if( ++i<argc ) params.kClust_cmd = std::string(argv[i]);
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
		
		if( arg=="-o" ){
			if( ++i<argc ){ 
				params.set_aln_db_dir( argv[i] );
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		
		if( arg=="-consdb" ){
			if( ++i<argc ){ 
				params.set_consensus_dbfile( argv[i] );
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		
		if( arg=="-cons" ){
			if( ++i<argc ){
				params.consensus_method = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		
		if( arg=="-p" ){
			params.parallel_aln_comp = true;
			continue;
		}


		if( arg=="--write-add-header" ){
			params.write_add_header = true;
			continue;
		}

		if( arg=="--no_advanced_clustering" ){
			params.advanced_clustering = false;
			continue;
		}
		
		if( arg=="--merge-ncbi-headers" ){
			params.consensus_header_merging_method = Clusters::NCBI_HEADER_MERGING;
			continue;
		}	
		if( arg=="--merge-uniprot-headers" ){
			params.consensus_header_merging_method = Clusters::UNIPROT_HEADER_MERGING;
			continue;
		}

		
		if( arg=="--no-pseudo-headers" ){
			params.write_pseudo_headers = false;
			continue;
		}

		if( arg=="-h" || arg=="-help" || arg=="--h" || arg=="--help" ){ 
			usage(std::cerr);
			exit(2);
		}
		throw MyException("Unknown command-line parameter found: '%s'", argv[i]);
	}
	if (params.get_aln_db_dir() == "") {
		params.set_aln_db_dir(params.get_working_dir() + "alignments");
	}
}

std::ostream& usage( std::ostream &out ){

	char buffer[512];
	sprintf(buffer, "Usage: ./kClust_mkAln -c '[command]' -d [directory] [options]\n");
	out << buffer;
	out << std::endl;

	out << "Version 1.0" << endl << endl;

	out << "kClust_mkAln creates alignment files and a consensus DB file." << endl;
	out << "Written by Christian Mayer (christian.eberhard.mayer@googlemail.com) and Maria Hauser (mhauser@genzentrum.lmu.de)." << endl;


	out << "Required arguments:" << endl;
	sprintf(buffer, " %-20s %-20s: %s\n", "-c", "'[command]'", "Command for multiple sequence alignment program (input/output must be in FastA format)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-d", "[directory]", "Directory with the result files of kClust.");
	out << buffer;
	out << std::endl;
	out << "Optional arguments:" << endl;
	sprintf(buffer, " %-20s %-20s: %s\n", "-cl", "[file]", "Location of the kClust binary for further clustering of the large clusters. (default=/cluster/toolkit/production/bioprogs/kClust/bin/kClust)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-p", "", "Parallelized calculation of the multiple sequence alignments.");
    out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-o", "[directory]", "Directory for the created alignments. (default: kClust-result-dir/alignments/)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-consdb", "[file]", "Filename for the consensus database (default: working-dir/consensus.fas)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-cons", "[0-5]", "Consensus calculation method. (default: 1)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "1: conventional method");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "2: maximum profile score with consensus");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "3: minimize the relative entropy between p(a) and P(a|x)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "4: minimize quadratic difference between profile score and matrix score");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "5: minimize expected score difference between homologs and random sequences");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--write-add-header", "", "Write consensus header line with # at the beginning of the FASTA alignment (default=none).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--merge-ncbi-headers", "", "Compress NCBI headers for consensus sequences (default=concatenation).");
	out << buffer;
	sprintf(buffer, " %-20s %-17s: %s\n", "--merge-uniprot-headers", "", "Compress Uniprot headers for consensus sequences (default=concatenation).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--no-pseudo-headers", "", "Write original header in cluster alignments (default=pseudo-headers).");
	out << buffer;
	out << std::endl;
	out << " Example: ./kClust_mkAln -c '/usr/bin/kalign -q -i $infile -o $outfile' -d /tmp/clustering/" << endl << endl;
	return out;
}
