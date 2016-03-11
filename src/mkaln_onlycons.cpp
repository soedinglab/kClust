/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 *                                                                         *
 *   rewritten by Michael Remmert (2008)                                   *
 *   remmert@lmb.uni-muenchen.de                                           *
 ***************************************************************************/
#include <boost/regex.hpp>

#include "simple_hash.h"
#include "alignment.h"
#include "mkaln.h"
#include "converter.h"

MkAln::MkAln(Params *params, Matrix *m) throw (std::exception):
	clu(params->get_nodes_dmp_file().c_str(), m),
	params(params){

	boost::regex expr("^(.*?)\\$infile(.*?)\\$outfile(.*)$");
	boost::cmatch matches;
	if( boost::regex_match(params->alignment_cmd.c_str(), matches, expr) ){
		cmd_first_part  = matches.str(1);
		cmd_second_part = matches.str(2);
		cmd_third_part  = matches.str(3);
	}else 
		throw MyException("Invalid alignment command '%s'", params->alignment_cmd.c_str());

}

MkAln::~MkAln(){}

void MkAln::make_alignments(	Fasta_db_reader *dbr, 
								const Clusters::header_merging_method_t merging_method
								) throw (std::exception){
	char buffer[10];
	buffer[9] = '\0';

	size_t error_alns = 0;
	
	std::ofstream consdb(params->get_consensus_dbfile().c_str());
		
	for( size_t i=start_index; i<=clu.dbsize; ++i){
		//skip members
		if( clu.representatives[i]==0 ) continue;
		++cl_idx;
		std::vector<Clusters::Accessions>  org_vec;
		Simple_hash<size_t>      org_check(1024);
		Simple_hash<Clusters::prio_string> key_check(1024);
		//used for simple concatenation
		Clusters::string_vector_type header_vec;

		Simple_list<Clusters::Member>::Iterator it  = clu.representatives[i]->begin();
		Simple_list<Clusters::Member>::Iterator end = clu.representatives[i]->end();

		Sequence *r         = dbr->get( (*it).idx, true, true);
		size_t cluster_size = clu.get_member_count( r->index() );
		Converter::get_file_name( buffer, r->index() );

		/* Comment out for only consensus-DB
		std::string cluster_filename     = params->get_working_dir() + std::string(buffer) + ".fas";
		std::string cluster_aln_filename = params->get_working_dir() + std::string(buffer) + ".afas";
		std::ofstream cluster_file;

		if( cluster_size>1 ){
			cluster_file.open(cluster_filename.c_str());
			r->replace_J();
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( cluster_file );
			} else {
				r->write( cluster_file );
			}
		}
		*/

		if( merging_method==Clusters::SIMPLE_CONCATENATION ){
			header_vec.push_back( std::string(r->get_header()) );
		}else{
			clu.eval_header( r->get_header(), org_vec, org_check, key_check );	
		}

		++it;
		while( it!=end ){
			Sequence *s = dbr->get( (*it).idx, true, true);

			/* Comment out for only consensus-DB
			s->replace_J();
			if (params->write_pseudo_headers) {
				s->write_with_pseudo_header( cluster_file );
			} else {
				s->write( cluster_file );
			}
			*/

			if( merging_method==Clusters::SIMPLE_CONCATENATION ){
				header_vec.push_back( std::string(s->get_header()) );
			}else{
				clu.eval_header( s->get_header(), org_vec, org_check, key_check );	
			}
			delete s;
			++it;
		}

		//create consensus header
		std::stringstream header_stream;
		if( merging_method==Clusters::SIMPLE_CONCATENATION ){
			clu.create_merged_header_simple(header_stream, header_vec, r, cluster_size);
		}else{
			clu.create_merged_header_advanced(	header_stream, 
															r, 
															cluster_size, 
															org_vec, 
															org_check, 
															key_check);	
		}

		
		//std::string alignment_filename = params->get_aln_db_dir() + std::string(buffer) + ".fas";
		std::string alignment_filename = params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".fas";
		/* Comment out for only consensus-DB
		std::ofstream alignment_file(alignment_filename.c_str());
		*/
		
		consdb << header_stream.str();  // consensus header
		/* Comment out for only consensus-DB
		if (params->write_add_header) {
			alignment_file << "#" << header_stream.str();  // consensus header
		}
		*/

		if( cluster_size>1 ){
			/* Comment out for only consensus-DB
			int ret = 0;
			std::string cmd = "";
			cluster_file.close();
			if (cluster_size < 8000) {
				cmd = cmd_first_part + cluster_filename + cmd_second_part + cluster_aln_filename + cmd_third_part;
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
			} else {    	// more than 8000 sequences in cluster - prefiltering
				// create tmp directory
				std::string tmp_dir = params->get_working_dir() + "temp_kClust";
				cmd = "mkdir " + tmp_dir;
				ret = system( cmd.c_str() );
				if (ret != 0)
					throw MyException("Cannot create '%s'!", tmp_dir.c_str());
				// cluster sequences to 90%
				std::string cmd = "/cluster/user/michael/HHblast/tools/kClust/bin/kClust -i " + cluster_filename + " -d " + tmp_dir + " -s 4.74";
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
				// create consensus file
				cmd = "/cluster/user/michael/HHblast/tools/kClust/bin/kClust_mkAln -c '" + params->alignment_cmd + "' -d " + tmp_dir;
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
				// align consensus file
				cmd = cmd_first_part + tmp_dir + "/consensus.fas" + cmd_second_part + cluster_aln_filename + cmd_third_part;
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
				// clean tmp-directory
				system( ("rm -r " + tmp_dir).c_str() );
			}
			
			if( ret==0 ){
				Alignment aln( cluster_aln_filename.c_str(), Alignment::fasta, clu.matrix);
				std::stringstream aln_stream;
				aln.write( aln_stream, params->alignment_output_format );
				consdb << aln.get_consensus_sequence(80, params->consensus_method);  // consensus sequence
				alignment_file << aln_stream.str();        // cluster alignment
				if( remove( cluster_filename.c_str() ) )
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl; 
				if( remove( cluster_aln_filename.c_str() ) ) 
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl;
			}else{
				std::cerr << "      Error: External program terminated with non-zero exit state. (CMD: " << cmd << ")" << std::endl;
				++error_alns;
				std::stringstream tmpstream;
				r->write_sequence( tmpstream );
				
				consdb << tmpstream.str();   // consensus sequence
				tmpstream.str("");
				if (params->write_pseudo_headers) {
					r->write_with_pseudo_header( tmpstream, 80);
				} else {	
					r->write( tmpstream, 80);
				}
				alignment_file << tmpstream.str();        // cluster alignment
				// Comment out for DEBUG
				//if( remove( cluster_filename.c_str() ) )
				//	std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl;
				//
			}
			*/
			////////////////////////////
			// For only consensus DB:
			Alignment aln( alignment_filename.c_str(), Alignment::fasta, clu.matrix);
			std::stringstream aln_stream;
			consdb << aln.get_consensus_sequence(80, params->consensus_method);  // consensus sequence
			////////////////////////////
		}else{
			std::cerr << "      Singleton " << buffer << std::endl;
			std::stringstream tmpstream;
			r->write_sequence( tmpstream );
			consdb << tmpstream.str();   // consensus sequence
			tmpstream.str("");
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( tmpstream, 80);
			} else {	
				r->write( tmpstream, 80);
			}
			/* Comment out for only consensus-DB
			alignment_file << tmpstream.str();        // cluster alignment
			*/
		}
		
		/* Comment out for only consensus-DB
		alignment_file.close();
		*/
		delete r;		
	}
	if(error_alns!=0)
		std::cerr << "   There were " << error_alns << " error alignments." << std::endl;
		
	consdb.close();
}
