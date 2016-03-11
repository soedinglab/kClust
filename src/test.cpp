/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer   *
 *   christian.eberhard.mayer@googlemail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <mysql++.h>
#include "custom.h"


#include "mysql_spec_structures.h"


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
//#include <boost/regex.hpp>

#include "matrix.h"
#include "fasta_db_reader.h"
#include "base_table.h"
#include "chronometer.h"
#include "kmer.h"
#include "recycle_table.h"
#include "kmeraln.h"
#include "clusters.h"
#include "alignment.h"

#include "packed_sequence.h"
#include "simple_hash.h"

#include "ident_clust.h"
#include "params.h"
//requires -msse
//#include <emmintrin.h>
using std::string;
using std::ifstream;

void read_sw(const char*, int &Lq, int &dbsize, float **ptr);
void read_kmer(const char*, int Lq, int dbsize, float **ptr);
/**
 * 
 * @param argc 
 * @param argv[] 
 * @return 
 */



int main(int argc, const char **argv){
	if( argc!=2 ){
		std::cerr << "Usage: ./test [fasta-file]" << std::endl; 
		return 1;
	}
	
	//const size_t l = atoi(argv[1]);

	Chronometer chron;
	/*Params p;
	p.set_dbfile( std::string(argv[1]) );
	Matrix m( Matrix::blosum62 );
	Ident_clust ic(&p);
	ic.cluster();
	chron.print_program_memory_usage( std::cerr );
	ic.write_results();
	chron.print_program_memory_usage( std::cerr );
	return 0;
}*/

	Matrix m( Matrix::static_blosum62 );
	m.print( std::cout );


	Kmer k(3,4.4, &m);
	Fasta_db_reader dbr( argv[1], &m );
	Base_table table( dbr.get_sequence_count(), 3, &m );

	Sequence **data = new Sequence*[dbr.get_sequence_count()+1];

	while(dbr.has_next()){
		Sequence *s = dbr.get_next(true, false, false);
		data[s->index()] = s;	
		table.add_representative( s );
		delete s;
		//s->write(std::cout);
		//k.create_idents_list(s);
		//k.print_kmer_seed_list();
	}

	dbr.reset();

	while(dbr.has_next()){
		Sequence *s = dbr.get_next(true, false, false);
		k.create_idents_list(s);
		table.count_idents( &k );
		Base_table::BasicMatchIterator it = table.begin_iter_all();
		Base_table::BasicMatchIterator end = table.end_iter_all();
		while( it!=end ){
			std::cout << s->index() << " <-> " << it.get_sequence_index() << " " << it.get_score();
		
			Sequence *t = data[it.get_sequence_index()];
			std::cout << "   -> " << k.count_idents( s, t, 3)<<std::endl;
			
			++it;
		}
		delete s;

	}




	return 0;

}
/*	m.print_p_back(std::cerr);
	Kmer kmer(6, 4.4, &m);
	std::string seq;
	std::cin >> seq;
	Sequence s(">none", seq.c_str(), &m, 1);
	for( int i=0; i<seq.length()-6+1; ++i){
		size_t k = kmer.kmer2index( &(s.get_sequence()[i]) );
		kmer.print_kmer(k, std::cout);
		std::cout  << " " << k << std::endl;
	}
	//Alignment aln(argv[1], Alignment::fasta, &m);
	//std::cout << aln.get_consensus_sequence(80) << std::endl;
	return 0;
	}
	/*double A_eff=0.0;
	for( size_t i=0; i<21; ++i){ 
		std::cout << m.get_p_back()[i] << std::endl;
		A_eff += m.get_p_back()[i]*m.get_p_back()[i];
	}
	std::cout << "a_eff: " <<1.0/ A_eff << std::endl;
	return 0;
	Fasta_db_reader dbr( argv[1], &m );
	Sequence *q = dbr.get_next();
	
//	Base_table table( 2 , 6, &m );
//	table.add_representative( s );	
	Kmer kmer( 6, 4.4f, &m);
	kmer.create_kmer_list_fast( q );
	//std::cout << kmer.get_list_length() << std::endl;
	//kmer.create_kmer_list( q, m.get_matrix() );
	//std::cout << kmer.get_list_length() << std::endl;

//	Kmer::kmer_score sc;
//	kmer.get_kmer_score_with( q, sc);
//	std::cerr << "Hash: score=" << sc.score << " count=" << sc.count << std::endl;
//	table.match( &kmer );
//	table.print_matches( std::cerr );
	std::cerr << "Memory db reader : " << dbr.get_memory_usage() << std::endl;
	std::cerr << "Memory sequence  : " << q->get_memory_usage() << std::endl;
	std::cerr << "Memory kmer      : " << kmer.get_memory_usage() << std::endl;
	size_t sum = dbr.get_memory_usage();
	sum += kmer.get_memory_usage();
	sum += q->get_memory_usage();
	std::cerr << "Sum              : " << sum/1024 << std::endl;
	chron.print_program_memory_usage(std::cerr);
	return 0;
	}
	//Kmeraln aln( 4, &m, 1.00e-3f, 50, 12.0f, 2.0f, 0.27456f );

	
/*
	int bins[10000];
	for( size_t i=0; i<10000; ++i ) bins[i]=0;
	int max_bin=0;

	Chronometer chron;
	chron.start();
	//m.round_bit_scores();
	for( size_t i=0; i<10; ++i ){

		//Sequence *query = new Sequence( l, &m);
		//query->write(std::cerr);
		//kmer.create_kmer_list_fast( query );

		//Sequence *templ = new Sequence( l, &m);
		//templ->write(std::cerr);
		table.reset();
		for( size_t l=0; l<10000; ++l ){
			kmer.create_kmer_list_fast( seqs[l] );
			table.match( &kmer ); 
		}
		/*table.add_representative( templ );
		table.match_full( &kmer );
		if( table.begin()!=table.end() ){
			aln.align( Kmeraln::FAST_ADDR, table.begin().get_match_list(), query, templ, m.get_matrix() );
			//std::cerr << aln.get_score() << std::endl;
			int at = (int)(aln.get_score()+0.5);
			if( at<0 ){
				std::cerr << "Score is negative! "<< std::endl;	
			}else{
				if( at>=10000 ){
					std::cerr << "Score is too high!!!" << std::endl;
				}else{
					++bins[at];
					if( at>max_bin ) max_bin=at;
			}
		}

		delete query;
		delete templ;
		}
	}
	chron.stop();
	std::cerr << "TIME: " << chron.getSeconds() << std::endl;
	
	for( size_t i=0; i<=max_bin; ++i ){
		std::cout << i << " " << bins[i] << " " << log(bins[i]+0.1) << std::endl;
	}
	return 0;
}
		
		//Fasta_db_reader dbreader( argv[2], &m );	
		//Base_table bt(dbreader.get_sequence_count(), 6, &m);	

		/*std::cerr << "Creating kmer-list..." << std::endl;
		Kmer kmer( 4, 3.0f, &m );
		for( size_t i=0; i<10; ++i )
			kmer.create_kmer_list( query, m.get_matrix() );
		//kmer.print_kmer_seed_list();
		std::cerr << "Number of kmers:" << kmer.get_list_length() << std::endl;
		/*
		std::cerr << "Reading database..." << std::endl;
		while( dbreader.has_next() ){
			Sequence *s = dbreader.get_next();
			bt.add_representative( s ); 
		}

		
		std::cerr << "Starting benchmark..." << std::endl;
		Chronometer c;
		c.start();
		for( size_t i=0; i<30; ++i ){
			bt.match( &kmer );
			size_t idx;
			float score;
			//bt.get_best_match( idx, score, query->length() );
			//std::cerr << "idx:" << idx << " score:" << score << std::endl;
		}
		c.stop();
		c.print_time( std::cerr );
		


		/*float *sw_scores = 0;
		int Lq           = 0;
		int dbsize       = 0;
		read_sw(argv[1], Lq, dbsize, &sw_scores);	
		float *kmer_scores;
		read_kmer(argv[2], Lq, dbsize, &kmer_scores);
		
		for( int i=0; i<dbsize; ++i){
			if( sw_scores[i]!=-1.0 )
				std::cout << sw_scores[i] << " " << kmer_scores[i] << std::endl;
		}*/


	/*union {
 		__m128 m128;
		float  f[4];
	} f1,f2,res;

	f1.m128   = _mm_set_ps(1.0,1.0,1.0,1.0);
	f2.m128   = _mm_set_ps(1.0,1.0,1.0,1.0);
	res.m128 = _mm_add_ps(f1.m128, f2.m128);
	std::cout << res.f[0] << " " << res.f[1] << " " << res.f[2] << " " << res.f[3] << std::endl;*/

/*
void read_kmer( const char *fn, int Lq, int dbsize, float **ptr ){
	ifstream in(fn);
	if( in.fail() ) throw MyException("Cannot open '%s'!\n", fn );
	float *scores = new float[dbsize];
	for( int i=0; i<dbsize; ++i) scores[i]=0.0f; 
	while( in.good() ){
		int index;
		float score;
		in >> index;
		in >> score;	
		scores[index-1] = score/Lq;
	}
	*ptr = scores;
}

void read_sw( const char *fn, int &Lq, int &dbsize, float **ptr ){
	ifstream in(fn);
	if( in.fail() ) throw MyException("Cannot open '%s'!\n", fn );
	char buffer[4096];
	while( in.good() ){
		if( in.get()!='#' ){
			in.unget();
			break;
		}else{
			in.getline(buffer, 4096);
			boost::regex expr("^(.+?):\\s+(\\d+).*$");
			boost::cmatch matches;

			if( regex_match(buffer, matches, expr) ){
				if( matches[1] == "LENGTH OF QUERY" )
					Lq = atoi(matches[2].first);
				else if( matches[1] == "SIZEOF DATABASE" ) 
					dbsize = atoi(matches[2].first);
			}
		}
	}
	std::cerr << "Length of query:" << Lq << std::endl;
	std::cerr << "Size of database:" << dbsize << std::endl;

	float *sw_scores = new float[dbsize];
	size_t cur=0;
	while( in.good() ){
		int index;
		int Lt;
		int aln_len;
		int equ_len; 
		float score;
		int idents;
		int pos_matches;
		int neg_zero_matches;
		int gaps;
		in >> index;
		in >> Lt;
		in >> aln_len;
		in >> equ_len;
		in >> score;
		in >> idents;
		in >> pos_matches;
		in >> neg_zero_matches;
		in >> gaps;

		int max_len = std::max(Lq, Lt);
		int min_len = std::min(Lq, Lt);
		double cov = min_len/(double)max_len;
		if( cov<0.8 )
			sw_scores[cur++] = -1.0f;
		else
			sw_scores[cur++] = score/Lq;
		//if( score/Lq > 3.0 ) 
		//std::cout << index << " : " << score/Lq << " " << score << " " << Lq << std::endl;
	}
	in.close();
	*ptr = sw_scores;
}
*/

