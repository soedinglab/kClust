/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <exception>


#include "k_clust.h"
#include "params.h"
#include "chronometer.h"
#include "my_exception.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::string;

ostream& usage( ostream &out );
void set_params( const int, const char**, Params& ) throw (std::exception);

int main(int argc, const char **argv){

	if( argc<5 ){
		usage(cerr);
		return 1;
	}

	try{
		Params params;
		set_params( argc, argv, params );
		params.check();

		K_Clust c( &params );

		if (params.profile_query){
			std::cerr << "Preprocessing...\n";
			c.prepare();
			std::cerr << "done.\n";
		}

		Chronometer chron;
		chron.start();

		c.cluster();

		//	c.kDP_benchmark_scop_tpfp();


		float cl_seconds = chron.getSnapShotSeconds();
		cerr << "Runtime for representatives: " << cl_seconds << " seconds" << endl;

		c.write_representatives_db();
		float wrep_seconds = chron.getSnapShotSeconds();
		cerr << "Runtime for write representatives: " << wrep_seconds - cl_seconds << " seconds" << endl;
		c.write_dmp_files();
		float wdmp_seconds = chron.getSnapShotSeconds();
		cerr << "Runtime for write dmp files: " << wdmp_seconds - wrep_seconds << " seconds" << endl;
//		c.write_debug_files((params.get_working_dir() +"debug.dat").c_str());
//		float wdebug_seconds = chron.getSnapShotSeconds();
//		cerr << "Runtime for write debug files: " << wdebug_seconds - wdmp_seconds << " seconds" << endl;

		if( params.refinement ){
			cerr << "Perform refinement" << endl;
			c.refine();
			cerr << "Runtime for refinement: " << chron.getSnapShotSeconds()-wdmp_seconds << " seconds" << endl;
		}
		//cerr << "After refinement - before free tables" << endl;
		c.free_tables_seedlists_kmeraln();

		if( params.refinement ){
			//cerr << "Write refined nodes" << endl;
			c.write_refined_nodes_dmp();
		}

		//cerr << "write debug files" << endl;
//		c.write_debug_files( (params.get_working_dir() +"debug_refined.dat").c_str() );

		//cerr << "cleanup tmp-files" << endl;
		c.cleanup_tmp_files();

		//cerr << "Print time" << endl;
		chron.stop();

		chron.print_time( std::cerr );

	}catch(const std::exception &e){
		std::cerr << std::endl << ">" << std::endl;
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << ">" << std::endl;
		return 3;
	}
	return 0;
}


void set_params(const int argc, const char **argv, Params &params) throw (std::exception){

	// Check -P option and set profile mode if necessary
	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-P" ){
			params.profile_query = true;
		}
	}

	// check for -s option and change corresponding values
	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-s" ){
			if( ++i<argc ){
				float s = atof(argv[i]);
				if (s == 0.0f ){
					s = 0.52f;
					params.clustering_threshold = 0.0f;
				}
				else{
					params.clustering_threshold = s;
				}


				// calculate sequence identity (page 58 thesis C.Mayer)
//				float seqid = (params.clustering_threshold + 0.68) / 0.06;
				float seqid = (s + 0.68) / 0.06;

				// adapt --filter-t option (page 23, 48 thesis C.Mayer)
				params.kmer_score_threshold = seqid * 0.02 - 0.05;

				// adapt --filter-T option (page 47 thesis C.Mayer)
				params.filter_kmer_similarity_threshold = seqid * 0.01 + 4.1;

				// adapt --kdp-T option
				params.kdp_kmer_similarity_threshold = seqid * 0.018 + 2.46;

				// prob of match hit
				params.kdp_p_m = exp(-2.193 * params.kdp_kmer_similarity_threshold + 11.88) / pow(21.0,4.0);

				// intra-diagonal gap penalty
				params.kdp_F = 6.85 * 2 * params.kdp_p_m / 0.05;

				if (params.profile_query)
					params.clustering_threshold = 0.0f;

			}else throw MyException("No value available for '%s'", arg.c_str());
		}
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-m4mer" ){
			if (params.profile_query){
				std::cerr << "\n>\nERROR: Usage of k-mer substitution matrices is not possible with -P option.\n>\n";
				exit(2);
			}
			if( ++i<argc ){
				string m4mer(argv[i]);
				params._4mer_m_file = m4mer;
			}else throw MyException("No value available for '%s'", arg.c_str());

			if (params._4mer_m_file != ""){
				params.filter_kmer_similarity_threshold -= 0.2f; //4.2f;
				params.kdp_kmer_similarity_threshold -= 0.5f; //2.5f;
				params.kmer_matrices = true;
			}
		}
		if( arg=="-m2mer" ){
			if (params.profile_query){
				std::cerr << "Usage of k-mer substitution matrices is not possible with -P option.\n";
				exit(2);
			}
			if( ++i<argc ){
				string m2mer(argv[i]);
				params._2mer_m_file = m2mer;
			}else throw MyException("No value available for '%s'", arg.c_str());

			if (params._2mer_m_file != ""){
				params.kmer_matrices = true;
			}
		}
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-sc" ){
			params.score_correction = true;
		}
	}


	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-m4mer" ){
			i++;
			continue;
		}
		if( arg=="-m2mer" ){
			i++;
			continue;
		}
		if( arg=="-sc" ){
			continue;
		}
		if( arg=="-i" ){
			if( ++i<argc ) params.set_input_fileordir( string(argv[i]) );

			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-P" ){
			// already set
			continue;
		}
		if( arg=="-R" ){
			params.refinement = true;
			continue;
		}
		if( arg=="-d" ){
			if( ++i<argc ) params.set_working_dir( string(argv[i]) );
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="-td" ){
			if( ++i<argc ) params.set_tmp_dir( string(argv[i]) );
			continue;
		}

		if( arg=="-M" ){
			if( ++i<argc ){
				long mb = atoi(argv[i]);
				params.memory_limit = 1024L*1024L*mb;
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="-s" ){
			if( ++i<argc ){
				//params.clustering_threshold = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="-e" ){
			if( ++i<argc ){
				params.clustering_evalue_threshold = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="-c" ){
			if( ++i<argc ){
				const float COV = atof(argv[i]);
				if( COV<=0.0 ){
					params.coverage_seq_len  = 0.0;
					params.aln_coverage_long = 0.0;
					params.check_coverage_seq_len = false;
				}else{
					params.coverage_seq_len  = COV;
					params.aln_coverage_long = COV;
					params.check_coverage_seq_len = true;
				}
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		/*   if( arg=="-c-short" ){
      if( ++i<argc ){ 
		params.aln_coverage_short = atof(argv[i]);
      }else throw MyException("No value available for '%s'", arg.c_str());
      continue;
    }*/

		if( arg=="--kdp-k" ){
			if( ++i<argc ){
				params.kdp_k = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--kdp-T" ){
			if( ++i<argc ){
				params.kdp_kmer_similarity_threshold = atof(argv[i]);
				// prob of match hit
				params.kdp_p_m = exp(-2.193 * params.kdp_kmer_similarity_threshold + 11.88) / pow(21.0,4.0);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}


		if( arg=="--kdp-F" ){
			if( ++i<argc ){
				params.kdp_F = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--kdp-G" ){
			if( ++i<argc ){
				params.kdp_G = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--kdp-E" ){
			if( ++i<argc ){
				params.kdp_E = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		/*
    if( arg=="--kdp-pm" ){
      if( ++i<argc ){ 
		params.kdp_p_m = atof(argv[i]);
      }else throw MyException("No value available for '%s'", arg.c_str());
      continue;
    }
		 */

		if( arg=="--kdp-delta" ){
			if( ++i<argc ){
				params.kdp_delta = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--filter-k" ){
			if( ++i<argc ){
				params.filter_k = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--filter-T" ){
			if( ++i<argc ){
				params.filter_kmer_similarity_threshold = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--filter-t" ){
			if( ++i<argc ){
				params.kmer_score_threshold = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}

		if( arg=="--merge-ncbi-headers" ){
			params.representatives_header_merging_method = Clusters::NCBI_HEADER_MERGING;
			continue;
		}

		if( arg=="--merge-uniprot-headers" ){
			params.representatives_header_merging_method = Clusters::UNIPROT_HEADER_MERGING;
			continue;
		}
		if (arg=="--write-time-benchmark"){
			params.write_time_bench = true;
			continue;
		}

		if( arg=="-h" || arg=="-help" || arg=="--h" || arg=="--help" ){
			usage(cerr);
			exit(2);
		}
		throw MyException("Unknown command-line parameter found: '%s'", argv[i]);
	}
	if (params.profile_query)
		params.clustering_threshold = 0.0f;
}

ostream& usage( ostream &out ){
	char buffer[512];
	sprintf(buffer, "Usage: ./kClust -i [fasta-db-file] -d [directory] [options]\n");
	out << buffer;
	out << endl;
	out << "Version 1.0" << endl << endl;

	out << "kClust is a clustering program for protein sequences." << endl;
	out << "Written by Christian Mayer (christian.eberhard.mayer@googlemail.com) and Maria Hauser (mhauser@genzentrum.lmu.de)" << endl << endl;

	out << "Required arguments:" << endl;

	sprintf(buffer, " %-20s %-20s: %s\n", "-i", "[fasta-db-file]", "Sequence database in fasta format or directory with the output of the previous kClust run if -P option is set.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-d", "[directory]", "Directory for temporary and result files.");
	out << buffer;

	out << std::endl << "Optional arguments: " << std::endl;
	sprintf(buffer, " %-20s %-20s: %s\n", "-M", "[megabytes]", "Memory limit for clustering (default=3500MB).");
	out << buffer;
//	sprintf(buffer, " %-20s %-20s: %s\n", "-R", " ", "Compute refined assignment of member sequences (slow!) (default=false).");
//	out << buffer;
//	sprintf(buffer, " %-20s %-20s: %s\n", "-m4mer", " ", "Location of the 4-mer substitution matrix for the calculation of the similar k-mer lists (default="", not using k-mer substitution matrices).");
//	out << buffer;
//	sprintf(buffer, " %-20s %-20s: %s\n", "-m2mer", " ", "Location of the 2-mer substitution matrix for the calculation of the similar k-mer lists (default="", not using k-mer substitution matrices).");
//	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-P", " ", "Cluster profiles computed from existing alignment files (default=false).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-sc", " ", "Use sequence background frequency score correction for the k-mer scores (default=false).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-td", "[directory]", "Directory for temporary files (default=WORKING_DIR/tmp)");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-s", "[float]", "Clustering threshold (score per column) (default=1.12 half bits ~ 30% sequence identity). Set to zero for the clustering based only on the e-value of the hit.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-e", "[float]", "Clustering E-value threshold (default=1.0e-4).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-c", "[float]", "Alignment coverage of the longer sequence (default=0.8).");
	out << buffer;
	//  sprintf(buffer, " %-20s %-20s: %s\n", "-c-short", "[float]", "Alignment coverage of the shorter sequence (default=0.8).");
	//  out << buffer;

	sprintf(buffer, " %-20s %-20s: %s\n", "--merge-ncbi-headers", " ", "Compress NCBI headers in representatives database, creating a merged header instead of the representative sequence header.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--merge-uniprot-headers", " ", "Compress Uniprot headers in representatives database, creating a merged header instead of the representative sequence header.");
	out << buffer;

	sprintf(buffer, " %-20s %-20s: %s\n", "--write-time-benchmark", " ", "Write time benchmark files, containing sequences which consume the most computation time (default=false).");
    out << buffer;

	out << std::endl << "Expert arguments: " << std::endl;

	sprintf(buffer, " %-20s %-20s: %s\n", "--filter-k", "[integer]", "Length of k-mers for similarity scoring filter (default=6).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--filter-T", "[float]", "Similarity threshold for filter k-mer generation (default=4.3 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--filter-t", "[float]", "k-mer score threshold for prefiltering (default=0.55 half bits).");
	out << buffer;

	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-k", "[integer]", "Length of k-mers for kDP alignments (default=4).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-T", "[float]", "Similarity threshold for kDP k-mer generation (default=2.9 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-G", "[float]", "Gap open penalty (default=12.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-E", "[float]", "Gap extension penalty (default=2.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-F", "[float]", "Intra-diagonal gap penalty (default=0.27 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-delta", "[integer]", "Width of delta window (default=50).");
	out << buffer;

	out << std::endl;
	out << "Sequence identity ~ score per column (see -s option):" << std::endl;
	out << "20%   30%   40%   50%   60%   70%   80%   90%   99%" << std::endl;
	out << "0.52  1.12  1.73  2.33  2.93  3.53  4.14  4.74  5.28" << std::endl << std::endl;

	out.flush();
	return out;
} 
