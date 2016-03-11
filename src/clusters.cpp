/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "clusters.h"


Clusters::Clusters( const size_t db_size, Matrix *m ):dbsize(db_size), matrix(m){
	representatives = new Simple_list<Member>*[dbsize+1];
	Simple_list<Member> **ptr = representatives;
	Simple_list<Member> **end = ptr+dbsize+1;
	while( ptr!=end ){ *ptr++=0; }
	belongs_to_block = new size_t[dbsize+1];
	_cluster_count = 0;
}

Clusters::Clusters( const char *clusters_dmp_file, Matrix *m ):matrix(m){
	std::cerr << "   Reading cluster dump file ... ";
	read_clusters_dmp_file(clusters_dmp_file);
	std::cerr << "ok" << std::endl;
}

const size_t Clusters::get_memory_usage()const{
	size_t ret = sizeof(Clusters);
	ret       += (dbsize+2)*sizeof(Simple_list<Member>*);
	ret       += (dbsize+2)*sizeof(Simple_list<Member>);
	ret       += (dbsize+2)*2*sizeof(size_t);
	return ret;
}

Clusters::~Clusters(){
	for(size_t i=0; i<=dbsize; ++i){
		delete representatives[i];
	}
	delete [] representatives;
	delete [] belongs_to_block;
}

void Clusters::create_cluster_for(const size_t index, const size_t blockindex ){
	representatives[index] = new Simple_list<Member>();
	representatives[index]->append( Member(index, 0.0f) );
	belongs_to_block[index] = blockindex;
	++_cluster_count;
}

size_t Clusters::get_number_of_clusters(){
	return _cluster_count;
}

void Clusters::add_member( const size_t idx_r, const size_t idx_m, const float score ){
	representatives[ idx_r ]->append( Member( idx_m, score) );
}

bool Clusters::is_representative_in( const size_t seqindex, const size_t blockindex){
	return ( (representatives[ seqindex ]!=0) && (belongs_to_block[seqindex]==blockindex) );
}

bool Clusters::is_representative( const size_t idx ){
	return (representatives[idx]!=0);
}

void Clusters::print(Fasta_db_reader *dbr){
	Sequence *s;
	for( size_t i=0; i<dbsize+1; ++i){
		if( representatives[i]!=0 ){
			Simple_list<Member>::Iterator it = representatives[i]->begin();
			Simple_list<Member>::Iterator end = representatives[i]->end();
			s = dbr->get( (*it).idx,false, false, true );
			std::cout << "Cluster: " << (*it).idx << " " <<  s->get_header() << std::endl;
			delete s;
			++it;
			while( it!=end ){
				s = dbr->get( (*it).idx, false, false, true );
				std::cout << "     Member: " << (*it).idx << " " << (*it).score << " " <<  s->get_header() <<  std::endl;
				delete s;
				++it;
			}
		}
	}
}

void Clusters::evaluate_members(){
	float *scores = new float[dbsize+1];
	float *ptr = scores;
	float *end = scores+dbsize+1;
	while( ptr!=end ) *ptr++=0.0f;
	size_t *best_rep = new size_t[dbsize+1];
	for( size_t i=1; i<dbsize+1; ++i ){
		if( representatives[i]!=0 ){
			Simple_list<Member>::Iterator it = representatives[i]->begin();
			while( it!=representatives[i]->end() ){
				float s = (*it).score;
				size_t idx = (*it).idx;
				if( s>scores[idx] ){
					scores[idx] = s;
					best_rep[idx] = i;
				}
				++it;
			}
		}
	}
	remove_members();
	//add members to best representative
	for( size_t i=1; i<dbsize+1; ++i ){
		if( representatives[i]==0 ){
			representatives[best_rep[i]]->append( Member(i, scores[i]) );
		}
	}
	delete [] scores;
	delete [] best_rep;
}

void Clusters::remove_members(){
	for( size_t i=0; i<=dbsize; ++i ){
		if( representatives[i]!=0 ){
			Simple_list<Member>::Iterator it = representatives[i]->begin();
			Member m( (*it).idx, 0.0f );
			delete representatives[i];
			representatives[i] = new Simple_list<Member>();
			representatives[i]->append( m );
		}
	}
}

void Clusters::write_dmp_files(	Fasta_db_reader *dbr, 
			const char *names, 
			const char *clusters,
			const bool profile_query) throw (std::exception){
	
	if (!profile_query) write_dmp_files_np(dbr, clusters, names);
	else write_dmp_files_profile(dbr, clusters, names);
}


void Clusters::write_dmp_files_np( Fasta_db_reader *dbr,
                        const char *clusters,
                        const char *names) throw (std::exception){

        std::ofstream names_dmp(names);
        if( !names_dmp ) throw MyException("Cannot open '%s' for writing!", names);

        std::ofstream clusters_dmp(clusters);
        if( !clusters_dmp ) throw MyException("Cannot open '%s' for writing!", clusters);

        clusters_dmp << "# " << dbsize << std::endl;

        for( size_t i=1; i<=dbsize; ++i){
                //skip members
                if( representatives[i]==0 ) continue;

                Simple_list<Member>::Iterator it  = representatives[i]->begin();
                Simple_list<Member>::Iterator end = representatives[i]->end();

                Sequence *r = dbr->get( (*it).idx, true, false);

                names_dmp    << r->index() << " " << r->get_header() << std::endl;
                clusters_dmp << r->index() << " " << r->index() << std::endl;
                ++it;

                while( it!=end ){
                        Sequence *s = dbr->get( (*it).idx, true, false);
                        names_dmp    << s->index() << " " << s->get_header() << std::endl;
                        clusters_dmp << s->index() << " " << r->index() << std::endl;
                        delete s;
                        ++it;
                }
                delete r;
        }
        names_dmp.close();
        clusters_dmp.close();
}

void Clusters::write_refined_nodes_dmp_file(const char *clusters) throw (std::exception){

        std::ofstream clusters_dmp(clusters);
        if( !clusters_dmp ) throw MyException("Cannot open '%s' for writing!", clusters);

        clusters_dmp << "# " << dbsize << std::endl;

        for( size_t i=1; i<=dbsize; ++i){
                //skip members
                if( representatives[i]==0 ) continue;

                Simple_list<Member>::Iterator it  = representatives[i]->begin();
                Simple_list<Member>::Iterator end = representatives[i]->end();
                size_t r_index = (*it).idx;
                clusters_dmp << r_index << " " << r_index << std::endl;
                ++it;
                while( it!=end ){
                        clusters_dmp << (*it).idx << " " << r_index << std::endl;
                        ++it;
                }
        }
        clusters_dmp.close();
}


void Clusters::write_dmp_files_profile(Fasta_db_reader *dbr, const char *clusters, const char *headers) throw (std::exception){
	int _singleton_count = 0;
	int _nsingleton_count = 0;
	std::ofstream clusters_dmp(clusters);
	if( !clusters_dmp ) throw MyException("Cannot open '%s' for writing!", clusters);

	std::ofstream headers_dmp;
	if (headers != NULL) {
		headers_dmp.open(headers);
		if( !headers_dmp ) throw MyException("Cannot open '%s' for writing!", headers);
	}
	
	clusters_dmp << "# " << dbr->get_orig_db_size() << std::endl;
	
	size_t * orig2idx = dbr->get_orig2idx();
	
	for( size_t j=1; j<=dbr->get_orig_db_size(); ++j){
		int count = 1;
		//skip members
//		std::cout << "orig_idx: " << j << ", idx: " << orig2idx[j] << "\n";
		size_t i = orig2idx[j];
		if( representatives[i]==0 ) continue;

		Simple_list<Member>::Iterator it  = representatives[i]->begin();
		Simple_list<Member>::Iterator end = representatives[i]->end();
		
		std::list<int> cluster_members;
		// add all the members of all the profiles in this cluster to the list
		while( it!=end ){
//			std::cout << "cluster member: " << (*it).idx << "\n";
			Profile* pp = dbr->get_profile((*it).idx, false, false, true);
			std::list<int> * members = pp->get_members();

			std::list<int>::iterator m_it = members->begin();
			std::list<int>::iterator m_end = members->end();
//			std::cout << "profile members:\n";
//			while(m_it != m_end){
//				std::cout << *m_it << " ";
//				m_it++;
//			}
//			std::cout << "\n";
			cluster_members.merge(*members);
			delete pp;
			++it;
		}
		
		int r = *(cluster_members.begin());
		
//		std::cout << "writing representative " << r << "\n" ;
		clusters_dmp << r << " " << r << "\n";
		if (headers != NULL) headers_dmp << r << " " << dbr->get_header_for(r) << "\n";
		
		cluster_members.remove(*(cluster_members.begin()));
		
		std::list<int>::iterator m_it = cluster_members.begin();
		std::list<int>::iterator m_end = cluster_members.end();
		
		while(m_it != m_end){
//			std::cout << "writing member " << *m_it << "\n" ;
			clusters_dmp << *m_it << " " << r << "\n";
			if (headers != NULL) headers_dmp << *m_it << " " << dbr->get_header_for(*m_it) << "\n";
			m_it++;
			count++;
		}
		if(count == 1) _singleton_count++;
		else _nsingleton_count++;
	}
	std::cout << "CLUSTERS: " << _cluster_count << "\n";
	std::cout << "SINGLETONS: " << _singleton_count << "\n";
	std::cout << "NOT SINGLETONS: " << _nsingleton_count << "\n";
	clusters_dmp.close();

}

void Clusters::write_representants_db( Fasta_db_reader *dbr, 
													const char *outfile, 
													const header_merging_method_t merging_method) 
													throw (std::exception){
	std::ofstream out(outfile);
	if( !out ) throw MyException("Cannot open '%s' for writing!", outfile);
	
	size_t * orig2idx = dbr->get_orig2idx();
	
	for( size_t j=1; j<=dbr->get_orig_db_size(); ++j){
		//skip members
		size_t i = orig2idx[j];
		//skip members
		if( representatives[i]==0 ) continue;
		
		//used for special header merging
		std::vector<Accessions> org_vec;
		Simple_hash<size_t> org_check(2048);
		Simple_hash<prio_string> key_check(2048);

		//used for simple concatenation
		string_vector_type header_vec;

		// iterator over the cluster members
		Simple_list<Member>::Iterator it  = representatives[i]->begin();
		Simple_list<Member>::Iterator end = representatives[i]->end();
		Sequence *r = dbr->get( (*it).idx, true, false, true );
		
		size_t cluster_size=0;

		if (merging_method!=Clusters::ORIGINAL_REPRESENTATIVE_HEADER){
			if( merging_method==Clusters::SIMPLE_CONCATENATION ){
				if (Profile* p = dynamic_cast<Profile*>(r)){
					for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
						header_vec.push_back( std::string(dbr->get_header_for(*it)));
				}
				else{
					header_vec.push_back( std::string(r->get_header()) );
				}
			}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
				if (Profile* p = dynamic_cast<Profile*>(r)){
					for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
						eval_ncbi_header( dbr->get_header_for(*it), org_vec, org_check, key_check );
				}
				else{
					eval_ncbi_header( r->get_header(), org_vec, org_check, key_check );
				}	
			}else{ // UNIPROT_HEADER_MERGING
				if (Profile* p = dynamic_cast<Profile*>(r)){
					for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
						eval_uniprot_header( dbr->get_header_for(*it), org_vec, org_check, key_check );
				}
				else{
					eval_uniprot_header( r->get_header(), org_vec, org_check, key_check );
				}
			}
			++it;
			cluster_size=1;

			while( it!=end ){

				// next cluster member
				Sequence *s = dbr->get( (*it).idx, true, false, true );

				if( merging_method==Clusters::SIMPLE_CONCATENATION ){
					if (Profile* p = dynamic_cast<Profile*>(s)){
						for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
							header_vec.push_back( std::string(dbr->get_header_for(*it)) );
					}
					else{
						header_vec.push_back( std::string(s->get_header()) );
					}
				}else if (merging_method==Clusters::NCBI_HEADER_MERGING){
					if (Profile* p = dynamic_cast<Profile*>(s)){
						for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
							eval_ncbi_header( dbr->get_header_for(*it), org_vec, org_check, key_check );
					}
					else{
						eval_ncbi_header( s->get_header(), org_vec, org_check, key_check );
					}
				}else{ // UNIPROT_HEADER_MERGING
					if (Profile* p = dynamic_cast<Profile*>(s)){
						for (std::list<int>::iterator it = (*p).get_members()->begin();  it != (*p).get_members()->end(); it++)
							eval_uniprot_header( dbr->get_header_for(*it), org_vec, org_check, key_check );
					}
					else{
						eval_uniprot_header( s->get_header(), org_vec, org_check, key_check );
					}
				}

				delete s;
				++cluster_size;
				++it;
			}
		}


		if (merging_method==Clusters::ORIGINAL_REPRESENTATIVE_HEADER){
			out << r->get_header() << "\n";
		} else{
			std::stringstream header_stream;
			if( merging_method==Clusters::SIMPLE_CONCATENATION ){
				create_merged_header_simple(header_stream, header_vec, r, cluster_size);
			}else{
				create_merged_header_advanced(header_stream, r, cluster_size, org_vec, org_check, key_check);
			}
			out << header_stream.str();
		}
		r->write_sequence( out );
		delete r;
	}
	out.close();
}

void Clusters::write_representative_numbers(const char *filename)throw (std::exception){
	std::ofstream out(filename);
	if( !out ) throw MyException("Cannot open '%s' for writing!", filename);
	for( size_t i=1; i<=dbsize; ++i){
		if( representatives[i]==0 ) continue;
		out << i << std::endl;
	}
	out.close();
}

void Clusters::write_debug_file(	Fasta_db_reader *dbr, const char *filename) throw (std::exception){

	std::ofstream out(filename);
	if( !out ) throw MyException("Cannot open '%s' for writing!", filename);

	for( size_t i=1; i<=dbsize; ++i){
		//skip members
		if( representatives[i]==0 ) continue;

		Simple_list<Member>::Iterator it  = representatives[i]->begin();
		Simple_list<Member>::Iterator end = representatives[i]->end();
		Sequence *r = dbr->get( (*it).idx, true, false, true );
		size_t memc = get_member_count( r->index() );
		out   << "Cluster: " << r->index() << " " << memc << " " << r->get_header() << std::endl;
		++it;

		while( it!=end ){
			Sequence *s = dbr->get( (*it).idx, true, false, true );
			out    << "   Member: " << s->index() << " score:" << (*it).score << " " << s->get_header() << std::endl;
			delete s;
			++it;
		}
		delete r;
	}
	out.close();
}

std::ostream& Clusters::create_merged_header_simple( std::ostream &hstream, string_vector_type &headers, Sequence* r, const size_t cluster_size )
{
    char buffer[10];
    buffer[9] = '\0';
    Converter::get_file_name( buffer, r->index() );
    hstream << ">cl|" << buffer << '|' << cluster_size << '|' << r->length() << " ";
    const char new_line_sep = 1;
    std::string tmp = headers[0];
    boost::trim(tmp);
    hstream << tmp;
    for( size_t i=1; i<std::min((int)headers.size(), 5); ++i ){
        tmp = headers[i];
        boost::trim(tmp);
        hstream << new_line_sep;
        hstream << tmp.substr(1);
    }
    hstream << std::endl;
    return hstream;
}

std::ostream& Clusters::create_merged_header_advanced( std::ostream &hstream, 
					      Sequence *r, 
					      const size_t cluster_size,
					      std::vector<Accessions> &org_vec,
					      Simple_hash<size_t> &org_check,
			      		Simple_hash<prio_string> &key_check ){

	const char new_line_sep = ' ';
	const char key_sep = '|';
	const char id_sep  = ',';

	std::vector<prio_string> key_vec;
	Simple_hash<prio_string>::Iterator hit = key_check.begin();
	while( hit!=key_check.end() ){
		key_vec.push_back( hit.getValue() );
		++hit;
	}
	
	//more weighting of keywords
	boost::regex expr_noninfo(".*?(hypothetical|unknown|putative|predicted|unnamed|probable|partial|possible|uncharacterized|fragment).*", boost::regex::icase);

	for( size_t j=0; j<key_vec.size(); ++j ){
		if( key_vec[j].prio==1 && boost::regex_match(key_vec[j].str.c_str(), expr_noninfo) ){
			key_vec[j].prio=0;
		}
	}

	std::sort(key_vec.begin(), key_vec.end(), prio_string::greater );
	
	char buffer[10];
	buffer[9] = '\0';
	Converter::get_file_name( buffer, r->index() );
	int max_prio=0;
	hstream << ">cl|" << buffer << '|' << cluster_size << '|' << r->length() << " ";
	for( size_t j=0; j<key_vec.size(); ++j ){
		if( key_vec[j].prio>max_prio ) max_prio = key_vec[j].prio;
		if( max_prio > key_vec[j].prio ) break;
		if (j>0) hstream << "; ";
		hstream << key_vec[j].str;
		
		//limit to 5 keywords	
		if( j==5 ) break;
	}	
	hstream << "." << new_line_sep;

	std::sort(org_vec.begin(), org_vec.end(), Accessions::greater);

	for( size_t j=0; j<org_vec.size(); ++j ){
		//print name of organism - do not remove the blanc or formatdb gets trouble
		hstream << "[" << org_vec[j].name << "]" << key_sep;
		//print '|gi,pdb,swissprot|'  accessions for the organism
		for( size_t k=0; k<org_vec[j].gis.size(); ++k ){
			if( org_vec[j].gis[k]!="" )
				hstream << org_vec[j].gis[k];
			if( org_vec[j].pdbs[k]!="" ) {
				if( org_vec[j].gis[k]!="" )
					hstream << id_sep;
				hstream <<  "pdb:" << org_vec[j].pdbs[k];
			}
			if( org_vec[j].sps[k]!="" )  {
				if( org_vec[j].gis[k]!="" || org_vec[j].pdbs[k]!="" )
					hstream << id_sep;
				hstream << "sp:" << org_vec[j].sps[k];
			}
			if (k < org_vec[j].gis.size()-1 ) 
				hstream << key_sep;
		}
//		if( j==200 ) break;
		//print pseudo newline if not processing last organism
		if( (j+1)!= org_vec.size() ) hstream << new_line_sep;
	}
	hstream << "." << std::endl;
	return hstream;
}

void Clusters::eval_ncbi_header( const char *header,
			    std::vector<Accessions> &org_vec,
			    Simple_hash<size_t> &org_check,
			    Simple_hash<prio_string> &key_check ) throw (std::exception){

	string_vector_type split_vec;

	// gi|51127382|emb|CAH17678.1| AV2 protein [Stachytarpheta leaf curl virus - [Hn31]]
		
	boost::regex capt_organism_expr("(.+?)\\[(.*(\\[.+\\])*.*)\\](.*?)");
	boost::regex expr_inc("^includes.*", boost::regex::icase);
	boost::regex expr_con("^contains.*", boost::regex::icase);
	boost::regex expr_sp("^second part.*", boost::regex::icase);
	boost::regex expr_fp("^first part.*", boost::regex::icase);

	boost::regex capt_gi_expr("^>*gi\\|(\\d+)(\\S+)\\s+(.+)$");
	boost::regex capt_pdb_expr("^.*?\\|pdb\\|(\\S{4})\\|(.|$).*");
	boost::regex capt_sp_expr("^.*?\\|sp\\|.+?\\|(.+?)(\\||$).*");

	//split headers by char \1
	boost::split( split_vec, header, boost::algorithm::is_from_range(1, 1) );
	for( size_t i=0; i<split_vec.size(); ++i ){
		int header_prio=1;
		boost::cmatch matches;
		if( boost::regex_match(split_vec[i].c_str(), matches, capt_gi_expr) ){

			std::string gi       = matches.str(1);
			std::string rest_acc = matches.str(2);
			std::string rest     = matches.str(3);

			boost::cmatch acc_matches;
			std::string pdb = "";
			std::string sp  = "";
			//extract pdb acc
			if( boost::regex_match(rest_acc.c_str(), acc_matches, capt_pdb_expr) ){
				header_prio = 3;
				pdb    = acc_matches.str(1);
				pdb   += acc_matches.str(2);
				//change chain identifier to lower case	
				if( pdb.length()>4 ) pdb[4] = tolower(pdb[4]);
				//std::cout << "  ->PDB match: " << pdb << std::endl;
			}
			//extract sp acc
			// swissprot headers have the greatest priority
			if( boost::regex_match(rest_acc.c_str(), acc_matches, capt_sp_expr) ){
				header_prio = 5;
				sp     = acc_matches.str(1);
			}

			boost::cmatch organism_matches;
			size_t index_of_org=0;
			if( boost::regex_match(rest.c_str(), organism_matches, capt_organism_expr)){

				std::string keys = organism_matches.str(1);
				keys            += organism_matches.str(3);

				std::string org  = organism_matches.str(2);
				if( 	!boost::regex_match(org.c_str(), expr_inc) &&
						!boost::regex_match(org.c_str(), expr_con) &&
						!boost::regex_match(org.c_str(), expr_fp)  &&
						!boost::regex_match(org.c_str(), expr_sp) ){

					std::string lower_case_org = org;
					boost::to_lower( lower_case_org );
					if( !org_check.exists( lower_case_org.c_str(), lower_case_org.length()) ){
						Accessions dummy;
						dummy.name = org;
						dummy.gis.push_back( gi );
						dummy.pdbs.push_back( pdb );
						dummy.sps.push_back( sp );
						dummy.prio = header_prio;
						org_check.put( lower_case_org.c_str(), lower_case_org.length(), org_vec.size() );
						org_vec.push_back( dummy );
					}else{
						index_of_org = *( org_check.get( lower_case_org.c_str(), lower_case_org.length() ) );
						org_vec[index_of_org].gis.push_back( gi );
						org_vec[index_of_org].pdbs.push_back( pdb );
						org_vec[index_of_org].sps.push_back( sp );
						org_vec[index_of_org].prio=header_prio;
					}	
				}else{
					//check if there are already gis with unknown organisms
					if( org_check.exists( "?", 1) ){
						index_of_org = *( org_check.get( "?", 1 ) );
						org_vec[index_of_org].gis.push_back( gi );
						org_vec[index_of_org].pdbs.push_back( pdb );
						org_vec[index_of_org].sps.push_back( sp );
						org_vec[index_of_org].prio=header_prio;
					}else{
						Accessions dummy;
						dummy.name = "?";
						dummy.gis.push_back( gi );
						dummy.pdbs.push_back( pdb );
						dummy.sps.push_back( sp );
						dummy.prio = header_prio;
						org_check.put( dummy.name.c_str(), dummy.name.length(), org_vec.size() );
						org_vec.push_back( dummy );
					}
				}
				
				boost::trim(keys);
				std::vector<std::string> header_split;
    			boost::split( header_split, keys, boost::is_any_of(";") );
    			
    			for( size_t k=0; k<header_split.size(); ++k ){
    				header_split[k] = header_split[k].substr(header_split[k].find("=")+1);
					std::string lower_case = header_split[k];
					boost::to_lower( lower_case );
					if( !key_check.exists( lower_case.c_str(), lower_case.length() ) ){ 
						key_check.put( lower_case.c_str(), lower_case.length(), prio_string(header_split[k], header_prio) );
					}
    			}

			}else{ //no organism in this header '[orgname]'
				if( org_check.exists( "?", 1) ){
						index_of_org = *( org_check.get( "?", 1 ) );
						org_vec[index_of_org].gis.push_back( gi );
						org_vec[index_of_org].pdbs.push_back( pdb );
						org_vec[index_of_org].sps.push_back( sp );
				}else{
					Accessions dummy;
					dummy.name = "?";
					dummy.gis.push_back( gi );
					dummy.pdbs.push_back( pdb );
					dummy.sps.push_back( sp );
					org_check.put( dummy.name.c_str(), dummy.name.length(), org_vec.size() );
					org_vec.push_back( dummy );
				}

				boost::trim(rest);
				std::vector<std::string> header_split;
    			boost::split( header_split, rest, boost::is_any_of(";") );
    			for( size_t k=0; k<header_split.size(); ++k ){
    				header_split[k] = header_split[k].substr(header_split[k].find("=")+1);
    				
					std::string lower_case = header_split[k];
					boost::to_lower( lower_case );
					if( !key_check.exists( lower_case.c_str(), lower_case.length() ) ){ 
						key_check.put( lower_case.c_str(), lower_case.length(), prio_string(header_split[k],header_prio) );
					}
    			}
			}
		}else{
			std::cerr << "WARNING! Header: " << header << " is not in NCBI-format!" << std::endl;
			// Throw exception only at small i values (problems in NR with large header lines)
			if (i < 100) { 
				throw MyException("Header '%s' is not in NCBI fasta-format!", header);
			}	
		}
	}
}

void Clusters::eval_uniprot_header( const char *header,
			    std::vector<Accessions> &org_vec,
			    Simple_hash<size_t> &org_check,
			    Simple_hash<prio_string> &key_check ) throw (std::exception){

	string_vector_type split_vec;

	boost::regex ncbi_header_expr("^>(\\S+)\\|(\\S+)\\|\\S+ (.+) OS=(.+?) (?:GN=.+ )?PE=\\d+ SV=\\d+$");
	int header_prio=1;
	boost::cmatch matches;
	if( boost::regex_match(header, matches, ncbi_header_expr) ){

		std::string db = matches.str(1);
		std::string id = matches.str(2);
		std::string descr = matches.str(3);
		std::string org = matches.str(4);
//		int seq_prio = atoi(matches.str(5).c_str()); // 1: best, 5: worst (= protein existence uncertain)

		boost::cmatch acc_matches;
		std::string sp = "";
		std::string tr = "";
		std::string pdb = "";
		//extract sp acc
		// swissprot headers have the greatest priority
		if( db.compare("sp") == 0 ){
			sp = id;
			header_prio = 3;
		}
		else{
			tr = id;
			header_prio = 1;
		}

		boost::cmatch organism_matches;
		size_t index_of_org = 0;

		std::string lower_case_org = org;
		boost::to_lower( lower_case_org );
		if( !org_check.exists( lower_case_org.c_str(), lower_case_org.length()) ){
			Accessions dummy;
			dummy.name = org;
			dummy.gis.push_back( tr );
			dummy.pdbs.push_back( pdb );
			dummy.sps.push_back( sp );
			dummy.prio = header_prio;
			org_check.put( lower_case_org.c_str(), lower_case_org.length(), org_vec.size() );
			org_vec.push_back( dummy );
		}else{
			index_of_org = *( org_check.get( lower_case_org.c_str(), lower_case_org.length() ) );
			org_vec[index_of_org].gis.push_back( tr );
			org_vec[index_of_org].pdbs.push_back( pdb );
			org_vec[index_of_org].sps.push_back( sp );
			org_vec[index_of_org].prio=header_prio;
		}

		boost::trim(descr);
		std::string lower_case = descr;
		boost::to_lower( lower_case );
		if( !key_check.exists( lower_case.c_str(), lower_case.length() ) ){
			key_check.put( lower_case.c_str(), lower_case.length(), prio_string(descr, header_prio) );
		}


	}else{
		throw MyException("Header '%s' is not in Uniprot fasta-format!", header);
	}
}

size_t Clusters::get_member_count( const size_t idx ){
	if( representatives[idx]!=0 ) return representatives[idx]->length();
	return 0;
}


std::ostream& Clusters::print_representants_with_members( std::ostream &out, Fasta_db_reader *dbr ){
	std::cerr << "Redundant sequences..." << std::endl;
	for( size_t i=1; i<=dbsize; ++i ){
		if( representatives[i]!=0 && representatives[i]->length()!=1 ){
			std::cerr << "r:" << i << " m:" << representatives[i]->length() << std::endl;
			Simple_list<Member>::Iterator it = representatives[i]->begin();
			while( it!=representatives[i]->end() ){
				size_t idx = (*it).idx;
				Sequence *s = dbr->get( idx, true, false, true );
				s->write( out );
				delete s;
				++it;
			}
		}
	}

	return out;
}



void Clusters::read_clusters_dmp_file( const char *filename ) throw (std::exception){
	std::ifstream in(filename);
	if(!filename) throw new MyException("Cannot open '%s'", filename);
	size_t node=0, parent=0;
	std::string comment;
	in >> comment;
	in >> dbsize;
	representatives = new Simple_list<Member>*[dbsize+1];
	Simple_list<Member> **ptr = representatives;
	Simple_list<Member> **end = ptr+dbsize+1;
	while( ptr!=end ){ *ptr++=0; }
	belongs_to_block = new size_t[dbsize+1];
	_cluster_count = 0;

	while( in.good() ){
		in >> node;
		if( in.eof() ) break;
		in >> parent;

		if( node==parent ){
			//representant
			create_cluster_for(node, 0);
		}else{
			//member
			add_member(parent, node, 0.0);
		}
	}	
	std::cout << " found " << (_cluster_count) << " clusters " << std::endl;
}

/*void Clusters::analyze_keywords(Fasta_db_reader *dbr) throw (std::exception){
	dbr->reset();
	Simple_hash<int> check(100000);

	boost::regex capt_organism_expr("(.+)\\[(.*(\\[.+\\])*.*)\\](.*?)");
	boost::regex capt_gi_expr("^>*gi\\|(\\d+)(\\S+)\\s+(.+)$");

	while( dbr->has_next() ){
		Sequence *s = dbr->get_next(false, false);
		string_vector_type split_vec;
		boost::split( split_vec, s->get_header(), boost::algorithm::is_from_range(1, 1) );
		for( size_t i=0; i<split_vec.size(); ++i ){
			boost::cmatch matches;
			if( boost::regex_match(split_vec[i].c_str(), matches, capt_gi_expr) ){

				std::string rest     = matches.str(3);
				boost::cmatch organism_matches;
				std::string lower_case_keys;
				if( boost::regex_match(rest.c_str(), organism_matches, capt_organism_expr)){
					std::string keys = organism_matches.str(1);
					keys            += organism_matches.str(3);
					boost::trim(keys);
					lower_case_keys = keys;
					boost::to_lower( lower_case_keys );
				}else{
					lower_case_keys = rest;
					boost::to_lower( lower_case_keys );
				}
				string_vector_type keys_split_vec;
				boost::split( keys_split_vec, lower_case_keys, boost::algorithm::is_space() );
				for( size_t j=0; j<keys_split_vec.size(); ++j){
					if( keys_split_vec[j].length()>2 ){
						if( !isalpha(keys_split_vec[j][0]) && !isdigit(keys_split_vec[j][0]) ){ 
							keys_split_vec[j] = keys_split_vec[j].substr(1);
						}
						const size_t len = keys_split_vec[j].length();
						if( !isalpha(keys_split_vec[j][len-1]) && !isdigit(keys_split_vec[j][len-1]) ){ 
							keys_split_vec[j] = keys_split_vec[j].substr(0, len-1);
						}


						if( check.exists(keys_split_vec[j].c_str(), keys_split_vec[j].length()) ){
							++(*check.get(keys_split_vec[j].c_str(), keys_split_vec[j].length()));
						}else{
							check.put(keys_split_vec[j].c_str(), keys_split_vec[j].length(), 1);
						}
					}
				}
			}
		}
		delete s;
	}

	Simple_hash<int>::Iterator it = check.begin();
	size_t keyword_count = 0;
	while( it!=check.end() ){
		if(it.getValue()>100) ++keyword_count;
		++it;
	}

	size_t p=0;
	prio_string *ar = new prio_string[keyword_count];
	it = check.begin();
	while( it!=check.end() ){
		if(it.getValue()>100){ 
			ar[p].str = (*it);
			ar[p].prio = it.getValue();
			++p;
		}
		++it;
	}

	std::sort(ar, ar+keyword_count);
	
	for( size_t i=0; i<keyword_count; ++i){
		std::cout << ar[i].str << " " << ar[i].prio << std::endl;
	}
	delete [] ar;
}*/
