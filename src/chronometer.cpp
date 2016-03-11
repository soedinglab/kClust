/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

#include "chronometer.h"

Chronometer::Chronometer(){ tstart = clock(); }

Chronometer::~Chronometer(){}

void Chronometer::start(){ tstart = clock(); }

void Chronometer::stop(){ tstop = clock(); }

long Chronometer::getClocks(){ return tstop-tstart; }

double Chronometer::getSeconds(){ return getClocks()/(double)CLOCKS_PER_SEC; }

double Chronometer::getSnapShotSeconds(){ return (clock()-tstart)/(double)CLOCKS_PER_SEC; }

std::ostream& Chronometer::print_time( std::ostream &out ){
	long t    = getClocks();
	long d    = t/(CLOCKS_PER_SEC*3600l*24l);
	t         = t%(CLOCKS_PER_SEC*3600l*24l);
	long h    = t/(CLOCKS_PER_SEC*3600l);
	t         = t%(CLOCKS_PER_SEC*3600l);
	long m    = t/(CLOCKS_PER_SEC*60l);
	t         = t%(CLOCKS_PER_SEC*60l);
	long s    = t/CLOCKS_PER_SEC;
	t         = t%CLOCKS_PER_SEC;
	double ms = 1000*(t/(double)CLOCKS_PER_SEC);
	out << "\nElapsed time:\n" << std::endl;
	char buffer[200];
	sprintf(buffer, "%lid %lih %lim %lis %4.2fms\n%-2.4f seconds\n", d, h, m, s, ms, getSeconds() );
	out << buffer;
	out.flush();
	return out;
}

std::ostream& Chronometer::print_program_memory_usage( std::ostream &out ){
	std::string fn = "/proc/self/status";
	std::ifstream in( fn.c_str() );
	if(!in){
		std::cerr << "Cannot open '" << fn << "'" << std::endl;
		return out;
	}
	out << "\nMemory usage:\n" << std::endl;
	std::string key;
	while(getline(in,key)){
		if( key.find("VmPeak")!=std::string::npos ){
			out<<key<<std::endl;
		}else if( key.find("VmSize")!=std::string::npos ){
			out<<key<<std::endl;
		}else if( key.find("VmHWM")!=std::string::npos ){
			out<<key<<std::endl;
		}else if( key.find("VmRSS")!=std::string::npos ){
			out<<key<<std::endl;
		}else if( key.find("VmData")!=std::string::npos ){
			out<<key<<std::endl;
		}
	}
	in.close();
	out.flush();
	return out;
}
