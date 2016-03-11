#ifndef CM_BLOCK_MANAGER_H
#define CM_BLOCK_MANAGER_H 1
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Maintains segments of a sequence database

#include <string>
#include <vector>

#include "simple_list.h"
#include "fasta_db_reader.h"
#include "clusters.h"

class BlockManager{

	public:
	
		// Block saves all the information about a block of representative sequences
		// i.e. block begin, end, table swap file and so on
		class Block{
			public:
				Block( const char *wdir=0, const int index=0, const size_t begin=0, const size_t end=0):
					begin(begin), end(end),index(index){
					char filename[200];
					if( wdir!=0 ) working_dir = std::string(wdir);
					sprintf(filename, "table%i.dat", index);
					file = std::string(filename);
				}

				bool operator==(const Block &b){
					return(b.index==index);
				}
			
				Block& operator=( const Block &rhs){
					begin = rhs.begin;
					end   = rhs.end;
					index = rhs.index;
					file  = rhs.file;
					length_of_shortest_seq = rhs.length_of_shortest_seq;
					length_of_longest_seq  = rhs.length_of_longest_seq;
					working_dir            = rhs.working_dir;
					return *this;
				}

				void set_index(const int index){ 
					this->index = index; 
					char filename[200];
					sprintf(filename, "table%i.dat", index);
					file = std::string(filename);
				}

				int get_index()const{ return index; }
				
				std::string get_filename(){ return working_dir+file; }

				std::ostream &print(std::ostream &out, const char* prefix)const{
					out << prefix << " segment: " << index                  << std::endl;
					out << prefix << "   begin: " << begin                  << std::endl;
					out << prefix << "     end: " << end                    << std::endl;
					out << prefix << "    file: " << working_dir+file       << std::endl;
					out << prefix << "   chars: " << characters             << std::endl;
					out << prefix << "     min: " << length_of_shortest_seq << std::endl;
					out << prefix << "     max: " << length_of_longest_seq  << std::endl << std::endl;
					return out;
				}
			public:
				size_t begin;
				size_t end;
				size_t characters;
				size_t length_of_shortest_seq;
				size_t length_of_longest_seq;
			private:
				int index;
				std::string file;
				std::string working_dir;
		};

		typedef Simple_list<Block>::Iterator Iterator;

		BlockManager(Fasta_db_reader *db_reader, Clusters *clusters);
		~BlockManager();

		Iterator begin();
		Iterator end();
		void free_sequences();
		void add_block( const Block &b );
		void save_representative( Sequence *r );
		
		Sequence* get_representative( size_t idx );
		void load( const Block &b );

		Block get_current_table_block();
		size_t get_current_table_index();

		size_t get_number_of_blocks();

		size_t get_memory_usage();
	
	private:
		const size_t dbsize;
		Fasta_db_reader *db_reader;
		Clusters        *clusters;
		Block           current_table_block;

		Simple_list<Block> block_list;

		size_t _current_sequence_memory_usage;
		Sequence **seqs;
};

#endif
