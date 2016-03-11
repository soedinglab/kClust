/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include <boost/regex.hpp>
#include "custom.h"
#include <mysql++.h>

#include "mysql_spec_structures.h"
#include "simple_hash.h"
#include "alignment.h"
#include "mysqldb.h"
#include "converter.h"

Mysqldb::Mysqldb(Params *params, Matrix *m) throw (std::exception):
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

	//set table name for corresponding structures
	mysql_sequence::table() = params->mysql_seq_table_name.c_str();
	mysql_cluster::table()  = params->mysql_clu_table_name.c_str();

	if( !params->continue_aln_comp ){
		start_index = 0;
		cl_idx      = 0;
		std::cerr << "   Init MySQL tables..." << std::endl;
		init_mysql_tables();
		std::cerr << "   ok" << std::endl;
	}else{
		recover();
	}
}

Mysqldb::~Mysqldb(){}

void Mysqldb::createdb(	Fasta_db_reader *dbr, 
								const Clusters::header_merging_method_t merging_method
								) throw (std::exception){
	char buffer[10];
	buffer[9] = '\0';

	size_t error_alns = 0;	
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

		Sequence *r         = dbr->get( (*it).idx, true, true, false );
		size_t cluster_size = clu.get_member_count( r->index() );
		Converter::get_file_name( buffer, r->index() );

		mysql_sequence seq_struct;
		seq_struct.id           = r->index();
		seq_struct.cluster      = cl_idx;
		seq_struct.header       = r->get_header();
		seq_struct.sequence     = r->get_sequence();
		insert_seq_into_db(seq_struct);

		std::string cluster_filename     = params->get_working_dir() + std::string(buffer) + ".fas";
		std::string cluster_aln_filename = params->get_working_dir() + std::string(buffer) + ".afas";
		std::ofstream cluster_file;

		if( cluster_size>1 ){
			cluster_file.open(cluster_filename.c_str());
			//r->replace_J();
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( cluster_file );
			} else {
				r->write( cluster_file );
			}
		}

		if( merging_method==Clusters::SIMPLE_CONCATENATION ){
			header_vec.push_back( std::string(r->get_header()) );
		}else{
			clu.eval_header( r->get_header(), org_vec, org_check, key_check );	
		}

		++it;
		while( it!=end ){
			Sequence *s = dbr->get( (*it).idx, true, true, false );

			seq_struct.id           = s->index();
			seq_struct.header       = s->get_header();
			seq_struct.sequence     = s->get_sequence();

			insert_seq_into_db(seq_struct);
			//s->replace_J();
			if (params->write_pseudo_headers) {
				s->write_with_pseudo_header( cluster_file );
			} else {
				s->write( cluster_file );
			}

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

		mysql_cluster cl;
		cl.id                 = cl_idx;
		cl.name               = buffer;
		cl.size               = cluster_size;
		cl.consensus_header = header_stream.str();

		if( cluster_size>1 ){
			cluster_file.close();
			std::string cmd = cmd_first_part + cluster_filename + cmd_second_part + cluster_aln_filename + cmd_third_part;
			std::cerr << "      Executing " << cmd << std::endl;
			int ret = system( cmd.c_str() );
			
			if( ret==0 ){
				Alignment aln( cluster_aln_filename.c_str(), Alignment::fasta, clu.matrix);
				std::stringstream aln_stream;
				aln.write( aln_stream, params->alignment_output_format );
				cl.consensus_sequence = aln.get_consensus_sequence(80);
				cl.alignment          = aln_stream.str();
				if( remove( cluster_filename.c_str() ) )
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl; 
				if( remove( cluster_aln_filename.c_str() ) ) 
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl;
			}else{
				std::cerr << "      Error: External program terminated with non-zero exit state. Cluster id:  " << cl.id << std::endl;
				++error_alns;
				std::stringstream tmpstream;
				r->write_sequence( tmpstream );
				cl.consensus_sequence = tmpstream.str();
				tmpstream.str("");
				if (params->write_pseudo_headers) {
					r->write_with_pseudo_header( tmpstream, 80);
				} else {	
					r->write( tmpstream, 80);
				}
				cl.alignment          = tmpstream.str();
				if( remove( cluster_filename.c_str() ) )
					std::cerr << "      Warning: Cannot remove " << cluster_filename << std::endl; 
			}
		}else{
			std::cerr << "      Singleton " << buffer << std::endl;
			std::stringstream tmpstream;
			r->write_sequence( tmpstream );
			cl.consensus_sequence = tmpstream.str();
			tmpstream.str("");
			if (params->write_pseudo_headers) {
				r->write_with_pseudo_header( tmpstream, 80);
			} else {	
				r->write( tmpstream, 80);
			}
			cl.alignment          = tmpstream.str();
		}
		insert_cluster_into_db(cl);
		delete r;		
	}
	if(error_alns!=0)
		std::cerr << "   There were " << error_alns << " error alignments." << std::endl;
}

void Mysqldb::insert_seq_into_db( mysql_sequence &mseq_struct) throw (std::exception){
	
	int retry=0;
	bool success = false;
	mysqlpp::Connection con;
	while(retry<4 && !success){
		++retry;
		con.connect(	params->mysql_db.c_str(), 
							params->mysql_host.c_str(), 
							params->mysql_user.c_str(), 
							params->mysql_passwd.c_str(), 0, 0, 3600);
		if( !con.connected() ){
			sleep(1000);
			continue;
		}
		while(retry<4 && !success){
			mysqlpp::Query  q = con.query();
			q.replace(mseq_struct);
			if( !q.execute().success )
				std::cerr << "Cannot execute MySQL statement, retrying..." << std::endl;
			else
				success = true;
		}
	}
	if( !success ) 
		throw MyException("Cannot execute MySQL statement for sequence %i", mseq_struct.id);
	con.close();
}

void Mysqldb::insert_cluster_into_db( mysql_cluster &cl ) throw (std::exception){
	int retry=0;
	bool success = false;
	mysqlpp::Connection con;
	while(retry<4 && !success){
		++retry;
		con.connect(	params->mysql_db.c_str(), 
							params->mysql_host.c_str(), 
							params->mysql_user.c_str(), 
							params->mysql_passwd.c_str(), 0, 0, 3600);
		if( !con.connected() ){
			sleep(1000);
			continue;
		}
		while(retry<4 && !success){
			mysqlpp::Query  q = con.query();
			q.replace(cl);
			if( !q.execute().success )
				std::cerr << "Cannot execute MySQL statement, retrying..." << std::endl;
			else
				success = true;
		}
	}
	if( !success ) 
		throw MyException("Cannot execute MySQL statement for cluster %i (%s)", cl.id, cl.name.c_str());
	con.close();
}

void Mysqldb::init_mysql_tables() throw (std::exception){	
	mysqlpp::Connection con(	params->mysql_db.c_str(), 
										params->mysql_host.c_str(), 
										params->mysql_user.c_str(), 
										params->mysql_passwd.c_str() );
	mysqlpp::Query  q = con.query();
	std::cerr << "      Dropping table '" << params->mysql_clu_table_name << "' if exists" << std::endl;
	if( !q.exec("DROP TABLE IF EXISTS " + params->mysql_clu_table_name) ) 
		throw MyException( "Cannot execute statement: " + q.str() );
	q.reset();
	q << "CREATE TABLE ";
	q << params->mysql_clu_table_name;
	q << " (id INT NOT NULL, ";
	q << "name VARCHAR(20), ";
	q << "size INT, ";
	q << "consensus_header LONGTEXT, "; //up to 4GB text
	q << "consensus_sequence LONGTEXT, "; //up to 4GB text
	q << "alignment LONGTEXT, "; //up to 4GB text
	q << "profil LONGTEXT, "; //up to 4GB text
	q << "PRIMARY KEY(id), UNIQUE KEY(name))";
	if( !q.execute().success ) throw MyException( "Cannot execute statement: " + q.str() );
	q.reset();
	
	std::cerr << "      Dropping table '" << params->mysql_seq_table_name << "' if exists" << std::endl;
	if( !q.exec("DROP TABLE IF EXISTS " + params->mysql_seq_table_name) ) 
		throw MyException( "Cannot execute statement: " + q.str() );
	q.reset();
	q << "CREATE TABLE ";
	q << params->mysql_seq_table_name;
	q << " (id INT NOT NULL, ";
	q << "cluster INT, ";
	q << "header LONGTEXT, "; //up to 4GB text
	q << "sequence LONGTEXT, "; //up to 4GB text
	q << "PRIMARY KEY(id), KEY(cluster))";
	if( !q.execute().success ) throw MyException( "Cannot execute statement: " + q.str() );
	con.close();
}

void Mysqldb::recover() throw (std::exception){
	mysqlpp::Connection con(	params->mysql_db.c_str(), 
										params->mysql_host.c_str(), 
										params->mysql_user.c_str(), 
										params->mysql_passwd.c_str());
	mysqlpp::Query  q = con.query();
	std::cerr << "      Checking cluster count..." << std::endl;
	q << "SELECT COUNT(*) FROM " + params->mysql_clu_table_name;
	mysqlpp::Result res = q.store();
 	const size_t cluster_count = ((unsigned int) res.at(0).at(0));
	std::cerr << "         Found " << cluster_count << " clusters." << std::endl;

	start_index = clu.get_number_of_clusters();
	size_t cls = 0;
	for( size_t i=0; i<=clu.dbsize; ++i){
		if( clu.representatives[i]==0 ) continue;
		++cls;
		if( cls==cluster_count ){
			start_index = i+1;
			cl_idx = cluster_count;
			break;
		}
	}
}

void Mysqldb::write_consensus_db_file() throw (std::exception){
	std::ofstream consdb(params->get_consensus_dbfile().c_str());
	if( !consdb )
		throw MyException( "Cannot open '%s' for writing!", params->get_consensus_dbfile().c_str() );
	mysqlpp::Connection con(	params->mysql_db.c_str(), 
										params->mysql_host.c_str(), 
										params->mysql_user.c_str(), 
										params->mysql_passwd.c_str());
	mysqlpp::Query  q = con.query();
	q << "SELECT COUNT(*) FROM ";
	q << params->mysql_clu_table_name;
	mysqlpp::Result fres = q.store();
 	const size_t cluster_count = ((unsigned int) fres.at(0).at(0));

	q.reset();
	q << "SELECT consensus_header,consensus_sequence FROM ";
	q << params->mysql_clu_table_name;

	mysqlpp::ResUse res = q.use();
	mysqlpp::Row row;
	for( size_t i=0; i<cluster_count; ++i){
		row = res.fetch_row();
		consdb << row.raw_data(0);
		consdb << row.raw_data(1);
	}
	consdb.close();
}

