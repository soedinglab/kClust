/***************************************************************************
 *   Copyright (C) 2007 by chris   *
 *   chris@Vega   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "sequence.h"

int main(int argc, char **argv){

	Matrix m(Matrix::static_blosum62);
	std::cout << "Matrix created." << std::endl;
	Sequence s(argv[1], Sequence::fasta, &m);
	std::cout << "Sequence created." << std::endl;
	m.print_int( std::cout );
	std::cout << std::endl << "Rescaled: " << std::endl;
	m.wn_rescale( &s );
	m.print_int( std::cout );
	std::cout << std::endl << "Error: " << m.get_error_sum() << std::endl;
	return 0;
}

