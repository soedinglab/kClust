#ifndef CM_MKALN_H
#define CM_MKALN_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 *                                                                         *
 *   rewritten by Michael Remmert (2008)                                   *
 *   remmert@lmb.uni-muenchen.de                                           *
 ***************************************************************************/

//DESCRIPTION
//Implements the setup of an alignment and consensus sequence database.

#include <stdio.h>
#include <vector>
#include <string>
#include "clusters.h"
#include "params.h"


class MkAln{
  
 public:
  MkAln(Params *params, Matrix *m) throw (std::exception);
  ~MkAln();
  void make_alignments(Fasta_db_reader *dbr, 
		       const Clusters::header_merging_method_t merging_method) throw (std::exception);
  void make_alignments_parallel(	Fasta_db_reader *dbr,
  								const Clusters::header_merging_method_t merging_method) throw (std::exception);
 private:
  Clusters clu;
  Params *params;
  std::string cmd_first_part;
  std::string cmd_second_part;
  std::string cmd_third_part;
  size_t start_index;
  size_t cl_idx;
};

#endif
