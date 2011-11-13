/*!
An example program for sfc++.

Copyright 2010, 2011 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "cstdlib"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "list"
#include "vector"

#include "sfc++.hpp"

using namespace std;
using namespace sfc;

int main(void)
{
	#ifndef SIZE
	#define SIZE 4
	#endif
	const unsigned int width = SIZE, height = SIZE, depth = SIZE;

	// test 2d version
	boost::array<unsigned int, 2> indices2d = boost::assign::list_of(width)(height);
	Sfc<2> mapping2d(indices2d);

	// precalculate sfc mapping of the whole grid
	mapping2d.cache_all();

	cout << "The cell at the lower left corner in the following will be the " << mapping2d.get_sfc_index(boost::assign::list_of(0)(0)) + 1
		<< "th cell visited by the space-filling curve."
		<< endl << endl;

	// test 2d index traversal
	cout << "Cell traversal order in a " << width << " by " << height << " grid:" << endl;
	unsigned int y = height - 1;
	while (true) {
		indices2d[1] = y;
		for (unsigned int x = 0; x < width; x++) {
			indices2d[0] = x;
			cout << setw(4) << setfill(' ') << mapping2d.get_sfc_index(indices2d);
		}
		cout << endl;

		if (y == 0) {
			break;
		} else {
			y--;
		}
	}
	cout << endl << endl;

	// test 2d sfc traversal
	ofstream traversal2d_file("2d.dat");
	for (unsigned int sfc_index = 0; sfc_index < mapping2d.size(); sfc_index++) {
		const boost::array<unsigned int, 2> indices = mapping2d.get_indices(sfc_index);
		traversal2d_file << indices[0] << "\t" << indices[1] << "\t" << sfc_index << endl;
	}
	traversal2d_file.close();

	cout << "2d traversal order was written into 2d.dat, visualize for example with gnuplot:" << endl
		<< "gnuplot -p -e \"set xrange [-1:" << SIZE
		<< "]; set yrange [-1:" << SIZE
		<< "]; plot '2d.dat' using 1:2 with lines\""
		<< endl << endl;


	// test 3d coordinate traversal
	boost::array<unsigned int, 3> indices3d = boost::assign::list_of(width)(height)(depth);
	Sfc<3> mapping3d(indices3d);

	// precalculate sfc mapping of the whole grid
	mapping3d.cache_all();

	ofstream traversal3d_file("3d.dat");
	for (unsigned int sfc_index = 0; sfc_index < mapping3d.size(); sfc_index++) {
		boost::array<unsigned int, 3> indices = mapping3d.get_indices(sfc_index);
		traversal3d_file << indices[0] << "\t"
			<< indices[1] << "\t"
			<< indices[2] << "\t"
			<< sfc_index
			<< endl;
	}
	traversal3d_file.close();

	cout << "3d traversal order was written into 3d.dat, visualize for example with gnuplot:" << endl
		<< "gnuplot -p -e \"set xrange [-1:" << SIZE
		<< "]; set yrange [-1:" << SIZE
		<< "]; set zrange [-1:" << SIZE
		<< "]; splot '3d.dat' using 1:2:3 with lines\""
		<< endl << endl;


	cout << "Change the grid size by defining SIZE when compiling, e.g. make SIZE=8" << endl;

	return EXIT_SUCCESS;
}
