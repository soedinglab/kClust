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
	start_index = 0;
}

MkAln::~MkAln(){}

void MkAln::make_alignments(	Fasta_db_reader *dbr, 
								const Clusters::header_merging_method_t merging_method
								) throw (std::exception){
	char buffer[10];
	buffer[9] = '\0';

	size_t error_alns = 0;
	
	std::ofstream consdb;
	if (params->consensus_method != 0){
		consdb.open(params->get_consensus_dbfile().c_str(), std::ios_base::app);
	}

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

		Sequence *r         = dbr->get( (*it).idx, true, false );
		size_t cluster_size = clu.get_member_count( r->index() );
		Converter::get_file_name( buffer, r->index() );
		
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

		if( merging_method==Clusters::SIMPLE_CONCATENATION ){
			header_vec.push_back( std::string(r->get_header()) );
		}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
			clu.eval_ncbi_header( r->get_header(), org_vec, org_check, key_check );
		}else{
			clu.eval_uniprot_header( r->get_header(), org_vec, org_check, key_check );
		}

		++it;
		int i = 1;
		while( it!=end ){
			Sequence *s = dbr->get( (*it).idx, true, false);

			s->replace_J();
			if (params->write_pseudo_headers) {
				s->write_with_pseudo_header( cluster_file );
			} else {
				s->write( cluster_file );
			}

			if( merging_method==Clusters::SIMPLE_CONCATENATION ){
				header_vec.push_back( std::string(s->get_header()) );
			}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
				clu.eval_ncbi_header( s->get_header(), org_vec, org_check, key_check );
			}else{
				clu.eval_uniprot_header( s->get_header(), org_vec, org_check, key_check );
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
		std::cout << "Filename: " << alignment_filename << "\n";
		std::ofstream alignment_file(alignment_filename.c_str());

		if (params->consensus_method != 0){
			consdb << header_stream.str();  // consensus header
		}
		if (params->write_add_header) {
			alignment_file << "#" << header_stream.str().substr(1);  // consensus header
		}

		if( cluster_size>1 ){
			int ret = 0;
			std::string cmd = "";
			cluster_file.close();
			if (cluster_size < 8000 || params->advanced_clustering == false) {
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
				std::string cmd = params->kClust_cmd + " -i " + cluster_filename + " -d " + tmp_dir + " -s 4.74";
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
				// for the case that there is only one cluster
				cmd = "cp " + tmp_dir + "/representatives.fas " + cluster_aln_filename;
				ret = system(cmd.c_str());
				// align representatives file
				cmd = cmd_first_part + tmp_dir + "/representatives.fas" + cmd_second_part + cluster_aln_filename + cmd_third_part;
				std::cerr << "      Executing " << cmd << std::endl;
				ret = system( cmd.c_str() );
				// clean tmp-directory
				int ret1 = system( ("rm -r " + tmp_dir).c_str() );
				if (ret1 != 0)
					throw MyException("Cannot remove '%s'!", tmp_dir.c_str());
			}

			if( ret==0 ){
				Alignment aln( cluster_aln_filename.c_str(), Alignment::fasta, clu.matrix);
				std::stringstream aln_stream;
				aln.write( aln_stream, params->alignment_output_format );
				if (params->consensus_method != 0){
					consdb << aln.get_consensus_sequence(80, params->consensus_method);  // consensus sequence
				}
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
				if( remove( cluster_filename.c_str() ) )
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl;

			}
		}else{
			std::cerr << "      Singleton " << buffer << std::endl;
			std::stringstream tmpstream;
			if (params->consensus_method != 0){
				r->write_sequence( tmpstream );
				consdb << tmpstream.str();   // consensus sequence
				tmpstream.str("");
			}
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( tmpstream, 80);
			} else {
				r->write( tmpstream, 80);
			}
			alignment_file << tmpstream.str();        // cluster alignment
		}

		alignment_file.close();
		delete r;
	}
	if(error_alns!=0)
		std::cerr << "   There were " << error_alns << " error alignments." << std::endl;

	consdb.close();
}

void MkAln::make_alignments_parallel(	Fasta_db_reader *dbr,
								const Clusters::header_merging_method_t merging_method
								) throw (std::exception){
	char buffer[10];
	buffer[9] = '\0';

	size_t error_alns = 0;

	std::ofstream consdb;
	if (params->consensus_method != 0){
		consdb.open(params->get_consensus_dbfile().c_str(), std::ios_base::app);
	}

	// generate cluster and header files
	std::cout << "Generating headers and temporary files...\n";
	for( size_t i=start_index; i<=clu.dbsize; ++i){
		if ((i*10)%clu.dbsize == 0){
			std::cout << ((i*100)/clu.dbsize) << "%... " << std::flush;
		}
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

		Sequence *r         = dbr->get( (*it).idx, true, false );
		size_t cluster_size = clu.get_member_count( r->index() );
		Converter::get_file_name( buffer, r->index() );

		// consensus header file
		std::string header_filename      = params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".hdr";
		std::ofstream cluster_file;
		std::ofstream header_file(header_filename.c_str());

		if (cluster_size == 1){
			// only one sequence in the cluster, no need to calculate an alignment
			// write the sequence in the end alignment file
			std::string cluster_filename     = params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".afas";
			cluster_file.open(cluster_filename.c_str());
			r->replace_J();
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( cluster_file );
			} else {
				r->write( cluster_file );
			}
		}
		else {
			// write sequences in a fasta file, it will be aligned in the next step
			std::string cluster_filename     = params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".clu";
			cluster_file.open(cluster_filename.c_str());
			r->replace_J();
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( cluster_file );
			} else {
				r->write( cluster_file );
			}
		}

		// add representative sequence header data
		if( merging_method==Clusters::SIMPLE_CONCATENATION ){
			header_vec.push_back( std::string(r->get_header()) );
		}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
			clu.eval_ncbi_header( r->get_header(), org_vec, org_check, key_check );
		}else{
			clu.eval_uniprot_header( r->get_header(), org_vec, org_check, key_check );
		}

		// add other sequences in the cluster
		++it;
		int i = 1;
		while( it!=end ){
			Sequence *s = dbr->get( (*it).idx, true, false);

			s->replace_J();
			if (params->write_pseudo_headers) {
				s->write_with_pseudo_header( cluster_file );
			} else {
				s->write( cluster_file );
			}

			// add header data
			if( merging_method==Clusters::SIMPLE_CONCATENATION ){
				header_vec.push_back( std::string(s->get_header()) );
			}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
				clu.eval_ncbi_header( s->get_header(), org_vec, org_check, key_check );
			}else{
				clu.eval_uniprot_header( s->get_header(), org_vec, org_check, key_check );
			}
			delete s;
			++it;
		}

		// calculate the consensus header
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

		if (params->consensus_method != 0){
			consdb << header_stream.str();  // consensus header
		}
		if (params->write_add_header) {
			header_file << "#" << header_stream.str().substr(1);  // consensus header
		}
		delete r;
		cluster_file.close();
		header_file.close();

		// cluster sequences in clusters with more than 8000 members and align only the representatives
		if (cluster_size > 8000){
			std::cout << "Advanced clustering of the cluster " << buffer << "!\n";
			std::string cluster_filename     = params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".clu";
			// create tmp directory
			std::string tmp_dir = params->get_working_dir() + "temp_kClust";
			std::string cmd = "mkdir " + tmp_dir;
			int ret = system( cmd.c_str() );
			if (ret != 0)
				throw MyException("Cannot create '%s'!", tmp_dir.c_str());

			// cluster sequences to 90%
			cmd = params->kClust_cmd + " -i " + cluster_filename + " -d " + tmp_dir + " -s 4.74";
			std::cerr << "      Executing " << cmd << std::endl;
			ret = system( cmd.c_str() );

			cmd = "cp " + tmp_dir + "/representatives.fas " + cluster_filename;
			ret = system(cmd.c_str());

			cmd = "cp " +cluster_filename  + " " + params->get_aln_db_dir() + buffer[0] + buffer[1] + "/" + std::string(buffer) + ".clu1";
			ret = system(cmd.c_str());

			// clean tmp-directory
			int ret1 = system( ("rm -r " + tmp_dir).c_str() );
			if (ret1 != 0)
				throw MyException("Cannot remove '%s'!", tmp_dir.c_str());
		}

	}
	consdb.close();
	std::cout << "100% done.\n Starting parallel alignments calculation...\n";

	// get the absolute path of the alignment directory
	std::string abs_cmd = "readlink -f " + params->get_aln_db_dir();
	FILE *fp = popen(abs_cmd.c_str(), "r");
	char buf[1024];
	fgets(buf, 1024, fp);
	fclose(fp);
	std::string aln_db_dir_abs_path(buf);
	aln_db_dir_abs_path = aln_db_dir_abs_path.substr(0, aln_db_dir_abs_path.length()-1);

	// calculate alignments
	std::string rsub_cmd = "rsub --logfile " + params->get_working_dir() + "kClust_mkAln_rsub.log -g '"
			+ aln_db_dir_abs_path + "/*/*.clu' --mult 1000  -c '"
			+ cmd_first_part + "FILENAME" + cmd_second_part + "DIRBASENAME.afas" + cmd_third_part + "'";
	std::cout << "Submitting alignment jobs:\n" + rsub_cmd + "\n";
	int ret = system( rsub_cmd.c_str() );
	if( ret!=0 ){
		std::cerr << "      Error: External program terminated with non-zero exit state. (CMD: " << rsub_cmd << ")" << std::endl;
	}

	// concatenate the consensus headers and the alignments into a single file
	std::string concat_cmd = "for i in " + aln_db_dir_abs_path + "/*/*.afas; do cat \"${i/.afas}\".hdr  $i > \"${i/.afas}\".fas; rm $i \"${i/.afas}\".hdr; done";
	std::cout << "Generating alignment files:\n" + concat_cmd + "\n";
	ret = system( concat_cmd.c_str() );
	if( ret!=0 ){
		std::cerr << "      Error: External program terminated with non-zero exit state. (CMD: " << concat_cmd << ")" << std::endl;
	}

	std::string rm_cmd = "for i in " + aln_db_dir_abs_path + "/*/*.clu; do rm $i; done";
	std::cout << "Removing temporary files:\n" + rm_cmd + "\n";
	ret = system( rm_cmd.c_str() );
	if( ret!=0 ){
		std::cerr << "      Error: External program terminated with non-zero exit state. (CMD: " << rm_cmd << ")" << std::endl;
	}
}
