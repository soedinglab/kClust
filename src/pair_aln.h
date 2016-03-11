#ifndef CM_ALN_RECONCILIATION_H
#define CM_ALN_RECONCILIATION_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION
//Container class for pairwise alignment information
//Implements residue-wise alignment comparison 

#include <fstream>
#include <iostream>
#include <vector>

class Pair_aln{

	public:
		Pair_aln(int);
		Pair_aln(const char *blastfile) throw (std::exception);
		~Pair_aln();
		size_t compare( Pair_aln& );
		size_t get_number_of_aligned_residues();
		std::ostream& print( std::ostream& );
		int *x;
		int *y;
		float score;
		size_t len_x;
		size_t aln_len_x;
		size_t len_y;
		size_t aln_len_y;
		size_t len_aln;
		size_t len_equ;
		size_t idents;
		size_t pos_matches;
		size_t neg_zero_matches;
		size_t gaps;
		size_t gaps_x;
		size_t gaps_y;
		size_t x_aln_start;
		size_t x_aln_end;
		size_t y_aln_start;
		size_t y_aln_end;

	private:
		void _eval_seqs(	std::vector<int> &xvec, 
								std::vector<int> &yvec, 
								int xs, 
								int ys, 
								const std::string xseq, 
								const std::string yseq );
};

#endif
