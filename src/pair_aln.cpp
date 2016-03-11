/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

#ifndef CM_NO_BOOST
#include <boost/regex.hpp>
#endif

#include "my_exception.h"
#include "pair_aln.h"

Pair_aln::Pair_aln(int len):
	score(0.0f),
	len_x(0),
	aln_len_x(0),
	len_y(0),
	aln_len_y(0),
	len_aln(0),
	len_equ(0),
	idents(0),
	pos_matches(0),
	neg_zero_matches(0),
	gaps(0),
	gaps_x(0),
	gaps_y(0),
	x_aln_start(0),
	x_aln_end(0),
	y_aln_start(0),
	y_aln_end(0){
	x = new int[len];
	y = new int[len];
}

Pair_aln::Pair_aln(const char *blastfile) throw (std::exception):
	score(0.0f),
	len_x(0),
	aln_len_x(0),
	len_y(0),
	aln_len_y(0),
	len_aln(0),
	len_equ(0),
	idents(0),
	pos_matches(0),
	neg_zero_matches(0),
	gaps(0),
	gaps_x(0),
	gaps_y(0),
	x_aln_start(0),
	x_aln_end(0),
	y_aln_start(0),
	y_aln_end(0){

	#ifndef CM_NO_BOOST
	std::ifstream in(blastfile);
	if(!in) throw MyException("Cannot open '%s'!", blastfile);
	
	boost::regex blanc("^\\s*$");
	boost::regex query_expr("^Query:\\s+(\\d+)\\s+(\\S+)\\s+.*");
	boost::regex subjt_expr("^Sbjct:\\s+(\\d+)\\s+(\\S+)\\s+.*");
	
	const size_t L = 4*1024*1024;
	char buffer[L];
	size_t consecutive_blanc_lines=0;

	std::string x_seq="", y_seq="";
	std::vector<int> xvec, yvec;
	int xs=-1, ys=-1;
	bool found_data = false;
	while( in.good() ){
		in.getline(buffer, L);
		//std::cout << buffer << std::endl;
		if( boost::regex_match(buffer, blanc) ){
			if( ys!=-1 && xs!=-1 ){
				//std::cout << x_seq << std::endl;
				//std::cout << y_seq << std::endl <<std::endl;
				_eval_seqs(xvec, yvec, xs, ys, x_seq, y_seq);
			}	
			++consecutive_blanc_lines;
			if( consecutive_blanc_lines>2 && found_data ) break;
			xs=-1;
			ys=-1;
		}else{
			consecutive_blanc_lines=0;
		}
		boost::cmatch matches;
		if( boost::regex_match(buffer, matches, query_expr) ){
			xs     = atoi( matches.str(1).c_str() );
			x_seq  = matches.str(2);
			found_data = true;
		}else
		if( boost::regex_match(buffer, matches, subjt_expr) ){
			ys     = atoi( matches.str(1).c_str() );
			y_seq  = matches.str(2);
			found_data = true;
		}
		
	}

	len_aln = xvec.size();
	x = new int[len_aln];
	y = new int[len_aln];
	for( size_t i=0; i<len_aln; ++i ){
		x[i] = xvec[i]-1;
		y[i] = yvec[i]-1;
	}
	#else
	throw MyException("Constructor not available - compile with libboost.");
	#endif
}


void Pair_aln::_eval_seqs(std::vector<int> &xvec, std::vector<int> &yvec, int xs, int ys, const std::string xseq, const std::string yseq){
	for( size_t i=0; i<xseq.size(); ++i ){
		if( xseq.at(i)=='-' ){
			xvec.push_back(0);
			yvec.push_back(ys++);
		}else if( yseq.at(i)=='-' ){
			xvec.push_back(xs++);
			yvec.push_back(0);
		}else{
			xvec.push_back(xs++);
			yvec.push_back(ys++);
		}
	}
}

Pair_aln::~Pair_aln(){
	delete [] x;
	delete [] y;
}

size_t Pair_aln::compare( Pair_aln &other ){
	size_t ret     = 0;
	size_t other_i = 0;
	size_t i       = 0;
	while( i<len_aln && other_i<other.len_aln){
		if( x[i]==-1 || y[i]==-1 ){++i; continue;}
		if( other.x[other_i]==-1 || other.y[other_i]==-1 ){ ++other_i; continue;}
		if( x[i]==other.x[other_i] && y[i]==other.y[other_i]  ){
			++ret;
			++i;
			++other_i;
		}else{
			if( x[i]<other.x[other_i] ){
				++i;
				continue;
			}
			if( x[i]>other.x[other_i] ){
				++other_i;
				continue;
			}
			if( x[i]==other.x[other_i] ){
				++other_i;++i;
				continue;
			}
		}
	}
	return ret;
}

size_t Pair_aln::get_number_of_aligned_residues(){
	size_t ret=0;
	for( size_t i=0; i<len_aln; ++i ) if( x[i]!=-1 && y[i]!=-1 ) ++ret;
	return ret;
}

std::ostream& Pair_aln::print(std::ostream &out){
	for ( size_t i=0; i<len_aln; ++i) out << i+1 << " " << x[i] << " " << y[i] << std::endl;
	return out;
}
