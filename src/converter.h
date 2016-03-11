#ifndef CM_CONVERTER_H
#define CM_CONVERTER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
//#include <algorithm>


class Converter{
	public:
		// writes the name of the alignment file into the buffer, for the representative sequence with index idx 
		// (just the name, not the whole path)
		static void get_file_name(char buffer[], size_t idx);
		// returns the index of the representative sequence of the cluster which alignment is saved in the file buffer
		static int get_idx_for_file(const char buffer[]);
	private:
		static char i2v [];
		static char i2k [];
		static int char2i [];
};

#endif
