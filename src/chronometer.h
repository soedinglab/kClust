#ifndef CM_CHRONOMETER_H
#define CM_CHRONOMETER_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//stopwatch

#include <ctime>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

class Chronometer{
	public:
		Chronometer();
		~Chronometer();
		long getClocks();
		double getSeconds();
		double getSnapShotSeconds();
		void start();
		void stop();
		unsigned long get_CLOCKS_PER_SEC(){return CLOCKS_PER_SEC;} 
		std::ostream& print_program_memory_usage( std::ostream &out );
		std::ostream& print_time( std::ostream &out );
		
	private:
		std::clock_t tstart;
		std::clock_t tstop;
};

#endif
