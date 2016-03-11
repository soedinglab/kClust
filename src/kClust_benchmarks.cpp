#include <iostream>
#include <cstdlib>
#include <fstream>
#include <exception>

#include "params.h"
#include "k_clust.h"
#include "my_exception.h"

std::ostream& usage( std::ostream &out );
void set_params( const int, const char**, Params& ) throw (std::exception);

int main(int argc, const char **argv){
	if( argc<5 ){
		usage(std::cerr);
		return 1;
	}

	try{
		Params params;
		set_params( argc, argv, params );

		K_Clust c( &params );
		if (params.benchmark == 1){
			c.pref_benchmark_scop_tpfp();
		}
		else if (params.benchmark == 2){
			c.pref_benchmark_scop_scatterplot();
		}
		else if (params.benchmark == 3){
			c.kDP_benchmark_scop_tpfp();
		}
		else if (params.benchmark == 4){
			c.kDP_benchmark_scop_TMscore();
		}
	}catch(const std::exception &e){
		std::cerr << std::endl << ">" << std::endl;
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << ">" << std::endl;
		return 3;
	}
	return 0;
}
void set_params(const int argc, const char **argv, Params &params) throw (std::exception){

	bool in = false;
	bool out = false;

	if (argc < 7){
		usage(std::cerr);
		exit(2);
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-b" ){
			if( ++i<argc ) {
				params.benchmark = atoi(argv[i]);
			}
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
	}

	if (params.benchmark < 0 || params.benchmark > 4){
		std::cerr << "Wrong or no value for the benchmark.\n";
		usage(std::cerr);
		exit(2);
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-m4mer" ){
			if( ++i<argc ){
				string m4mer(argv[i]);
				params._4mer_m_file = m4mer;
			}else throw MyException("No value available for '%s'", arg.c_str());

			if (params._4mer_m_file != ""){
				params.filter_kmer_similarity_threshold = 4.2f;
				params.kdp_kmer_similarity_threshold = 2.5f;
				params.kmer_matrices = true;
			}
		}
		if( arg=="-m2mer" ){
			if( ++i<argc ){
				string m2mer(argv[i]);
				params._2mer_m_file = m2mer;
			}else throw MyException("No value available for '%s'", arg.c_str());

			if (params._2mer_m_file != ""){
				params.filter_kmer_similarity_threshold = 4.2f;
				params.kdp_kmer_similarity_threshold = 2.5f;
				params.kmer_matrices = true;
			}
		}
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-sc" ){
			params.filter_kmer_similarity_threshold = 4.3f;
			params.kdp_kmer_similarity_threshold = 2.8f;
			params.score_correction = true;
		}
	}

	for(int i=1; i<argc; ++i){
		string arg(argv[i]);
		if( arg=="-b" ){
			++i;
			continue;
		}
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
			if( ++i<argc ) {
				if (params.benchmark == 4){
					params.input_dir = (string(argv[i]));
				}
				else{
					params.set_input_fileordir(string(argv[i]));
				}
				in = true;
			}
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-d" ){
			if( ++i<argc ) {
				params.set_working_dir( string(argv[i]) );
				out = true;
			}
			else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-id" ){
			if (params.benchmark == 1 || params.benchmark == 3){
				params.count_idents = true;
			}
			else {
				std::cerr << "Counting identities only possible with benchmark methods 1 and 3.\n";
				usage(std::cerr);
				exit(2);
			}
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
		if( arg=="--filter-T" ){
			if( ++i<argc ){
				params.filter_kmer_similarity_threshold = atof(argv[i]);
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
		if( arg=="-h" || arg=="-help" || arg=="--h" || arg=="--help" ){
			usage(std::cerr);
			exit(2);
		}
		throw MyException("Unknown command-line parameter found: '%s'", argv[i]);
	}
	if (!in || !out){
		usage(std::cerr);
		exit(2);
	}
	params.check_coverage_seq_len = false;
}

std::ostream& usage( std::ostream &out ){
	char buffer[512];
	sprintf(buffer, "Usage: ./kClust_benchmarks -b [benchmark-method] -i [fasta-db-file] -d [directory] [options]\n");
	out << buffer;
	out << "\n";

	out << "Required arguments:\n";

	sprintf(buffer, " %-20s %-20s: %s\n", "-b", "[1-4]", "Benchmark method:");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "1: Write prefiltering scores for all against all comparison on a SCOP database.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "2: Write prefiltering scores only for homologous sequences from the SCOP database (homologous = in the same superfamily).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "3: Write kDP scores for all against all comparison on a SCOP database.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "", "", "4: Write alignments for all-against-all comparison within given SCOP families (family files have to end with \"_db\").");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-i", "[fasta-db-file]", "SCOP database in fasta format or directory with SCOP families.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-d", "[directory]", "Directory for result files.");
	out << buffer;


	out << std::endl << "Optional arguments: " << std::endl;
	sprintf(buffer, " %-20s %-20s: %s\n", "-m4mer", " ", "Location of the 4-mer substitution matrix for the calculation of the similar k-mer lists (default="", not using k-mer substitution matrices).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-m2mer", " ", "Location of the 2-mer substitution matrix for the calculation of the similar k-mer lists (default="", not using k-mer substitution matrices).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-sc", " ", "Use sequence background frequency score correction for the k-mer scores (default=false).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-id", " ", "Count identities instead of similar k-mers. Only available for the benchmark method 1 (k=6).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-F", "[float]", "Intra-diagonal gap penalty (default=0.04 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-G", "[float]", "Gap open penalty (default=12.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-E", "[float]", "Gap extension penalty (default=2.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--filter-T", "[float]", "Similarity threshold for filter k-mer generation (default=4.4, -m=2.5, -sc=2.8 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--kdp-T", "[float]", "Similarity threshold for kDP k-mer generation (default=3.0, -m=2.5, -sc=2.8 half bits).");
	out << buffer;
	out << std::endl;
	out.flush();
	return out;
} 
