/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

#include <vector>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include "fasta_db_reader.h"
#include "recycle_table.h"
#include "kmeraln.h"

using std::endl;

std::string fasta1="";
std::string fasta2="";
float thr    = 2.9f;
float cov    = 0.0;
int k        = 4;
int delta    = 50;
float pm     = 1.00e-3; 
float G      = 12.0f; 
float E      = 2.0f; 
float F      = 0.27f;
//bool rescale = false;
bool roundm  = false;
bool scores  = false;
bool spaced  = false;
bool verbose = false;
bool exists(const char*);
std::ostream& usage(std::ostream &out);
void set_params( const int, const char** ) throw (std::exception);

int main(int argc, const char **argv){

	if( argc<2 ){
		usage(std::cerr);
		return 1;
	}
	
	try{
		set_params(argc, argv);
		if( spaced && k!=4 ) throw MyException("Spaced k-mers implemented for k=4 only.");

		Matrix m(Matrix::static_blosum62);
		if( roundm ) m.round_bit_scores();
	
		if( fasta2=="" ) fasta2 = fasta1;
		Fasta_db_reader *dbr = new Fasta_db_reader(fasta2.c_str(), &m, false, false);
	
		const size_t dbseqs_count = dbr->get_sequence_count()+1;
		Sequence **dbseqs = new Sequence*[dbseqs_count];
		Recycle_table table(dbr->get_sequence_count(), k, &m);
	
		while( dbr->has_next() ){
			Sequence *s = dbr->get_next(false, false);
			dbseqs[s->index()] = s;
			if( spaced )
				table.add_representative_spaced(s);
			else
				table.add_representative(s);
		}
		
		delete dbr;
		dbr = new Fasta_db_reader(fasta1.c_str(), &m, false, false);
	
		Kmer kmer(k, thr, false);
		Kmeraln kmeraln( k, &m, pm, delta, G, E, F );
		if( cov>0.0 ) kmeraln.set_coverage_criterion(cov);
		
		if(scores) std::cerr << "[sequence index] [sequence index] [alignment score] [score per column]" << std::endl;
		
		// TODO überprüfen, ob die Methode immer noch dasselbe macht nach dem rescale
		// TODO create_kmer_list Funktionen funktionieren so nur für Sequenzen
		while( dbr->has_next() ){
				Sequence *q =dbr->get_next(false, false);
				if( verbose ) std::cerr << "Processing " << q->index() << std::endl;
/*				if( rescale ){
					m.wn_rescale( q );
					if(spaced)
						kmer.create_kmer_list_spaced( q );
					else
						kmer.create_kmer_list( q );
				}else{*/
					if(spaced)  
						kmer.create_kmer_list_fast_spaced( q );
					else if( k==4 || k==6 ) 	
						kmer.create_kmer_list_fast( q );
					else
						kmer.create_kmer_list( q );
//				}
				table.match_full( &kmer );
				Recycle_table::MatchIterator        it = table.begin();
				const Recycle_table::MatchIterator end = table.end();
				while( it!=end ){
					if( cov>0.0){
						double coverage;
						if(q->length()>dbseqs[it.get_sequence_index()]->length() )
							coverage = dbseqs[it.get_sequence_index()]->length()/(double)q->length(); 
						else
							coverage = q->length()/(double)dbseqs[it.get_sequence_index()]->length(); 
					 	if(coverage < (cov-0.0001)){
							++it;
							continue;
						}		
					}
					kmeraln.align(Kmeraln::FAST_ADDR, it.get_match_list(), q, dbseqs[it.get_sequence_index()]);
					if(!scores) kmeraln.print(std::cout, q, dbseqs[it.get_sequence_index()], m.get_matrix() );
					else std::cout << q->index() << " " << it.get_sequence_index() << " " << kmeraln.get_score() << " " << kmeraln.get_rel_score() << std::endl;
					++it;	
				}
				delete q;	
		}
		
		for (size_t i=1; i<dbseqs_count; ++i) delete dbseqs[i];
		delete [] dbseqs;
		delete dbr;	
		
	}catch( std::exception &e){
		std::cerr << std::endl << ">" << std::endl;
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << ">" << std::endl;
		return 2;
	}
	return 0;
}


void set_params(const int argc, const char **argv) throw (std::exception){
	
	if( !exists(argv[1]) ) throw MyException("fasta-file-1: '%s' does not exist" , argv[1]);
	fasta1 = std::string(argv[1]);
 

	for(int i=2; i<argc; ++i){
		std::string arg(argv[i]);
		if( arg=="-c" ){
			if( ++i<argc ){ 
				cov = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-T" ){
			if( ++i<argc ){ 
				thr = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-k" ){
			if( ++i<argc ){ 
				k = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-D" ){
			if( ++i<argc ){ 
				delta = atoi(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-p" ){
			if( ++i<argc ){ 
				pm = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-G" ){
			if( ++i<argc ){ 
				G = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-E" ){
			if( ++i<argc ){ 
				E = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
		if( arg=="-F" ){
			if( ++i<argc ){ 
				F = atof(argv[i]);
			}else throw MyException("No value available for '%s'", arg.c_str());
			continue;
		}
/*		if( arg=="--rescale" ){
			rescale = true;
			continue;
		}*/
		if( arg=="--round" ){
			roundm = true;
			continue;
		}
		if( arg=="--scores" ){
			scores = true;
			continue;
		}
		if( arg=="--spaced" ){
			spaced = true;
			continue;
		}
		if( arg=="-v" ){
			verbose = true;
			continue;
		}

		if( arg=="-h" || arg=="-help" || arg=="--h" || arg=="--help" ){ 
			usage(std::cerr);
			exit(2);
		}

		if( fasta2=="" ){
			if( !exists(arg.c_str()) ) throw MyException("fasta-file-2: '%s' does not exist" , arg.c_str());
			fasta2 = arg;
			continue;
		}
		throw MyException("Unknown command-line parameter found: '%s'", argv[i]);
	}
}


std::ostream& usage( std::ostream &out ){
	char buffer[512];
	sprintf(buffer, "Usage: ./kDP [fasta-file-1] [[fasta-file-2] [options]]\n");
	out << buffer;
	out << endl;
	out << "Version 1.0" << endl << endl;

	out << "kDP is a program for all-against-all protein sequence comparison." << endl;

	out << "Optional arguments: " << endl;

	sprintf(buffer, " %-20s %-20s: %s\n", "[fasta-file-2]", "", "Compare all sequences of [fasta-file-1] to all sequences of [fasta-file-2]");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-T", "[float]", "Threshold for similar k-mers (default=2.9 half bits).");
	out << buffer;
//	sprintf(buffer, " %-20s %-20s: %s\n", "-k", "[int]", "Length of k-mers (default=4).");
//	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-c", "[float]", "Coverage (default=disabled).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-D", "[int]", "Width of delta-window (default=50).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-p", "[float]", "Probability for chance match (default=1e-3).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-G", "[float]", "Gap open penalty (default=12.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-E", "[float]", "Gap extension penalty (default=2.0 half bits).");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "-F", "[float]", "Gap penalty for intra-diagonal distance (default=0.27 half bits).");
	out << buffer;
//	sprintf(buffer, " %-20s %-20s: %s\n", "--rescale", "", "Rescale BLOSUM matrix for query amino acid frequencies with conjugate gradient method.");
//	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--round", "", "Round BLOSUM substitution scores to integers.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--scores", "", "Print scores only.");
	out << buffer;
	sprintf(buffer, " %-20s %-20s: %s\n", "--spaced", "", "Use spaced k-mers (XX0XX).");
	out << buffer;
	out << endl;
	out.flush();
	return out;
} 

bool exists(const char *fn){
	bool ret=false;
	struct stat fstats;
	if( stat(fn, &fstats )==0 )
		ret=true;
	return ret;
}
