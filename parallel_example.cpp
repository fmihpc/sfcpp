/*!
An example program for using sfc++ in parallel.

Copyright 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


For example on a AMD Phenom(tm) II X6 1075T this command:

make SIZE=300 && ./parallel_example

prints:

Serial mapping calculated in 69 seconds 
Parallel mapping calculated in 23 seconds using 6 thread(s)
*/

#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "list"
#include "utility"
#include "vector"

#include "sfc++.hpp"

#ifdef _OPENMP
#include "omp.h"
#else
int omp_get_num_threads(void)
{
	return 1;
}
#endif

using namespace std;
using namespace sfc;

int main(void)
{
	#ifndef SIZE
	#define SIZE 4
	#endif
	const unsigned int width = SIZE, height = SIZE, depth = SIZE;
	boost::array<unsigned int, 3> indices = boost::assign::list_of(width)(height)(depth);

	// cache grid <-> sfc mapping serially
	Sfc<3> mapping(indices);

	const time_t before_serial = time(NULL);
	mapping.cache_all();
	const time_t after_serial = time(NULL);

	cout << "Serial mapping calculated in " << after_serial - before_serial << " seconds " << endl;

	// write out sfc traversal order calculated serially for verification
	/*ofstream serial_file("serial.dat");
	for (unsigned int sfc_index = 0; sfc_index < mapping.size(); sfc_index++) {

		if (!mapping.is_cached(sfc_index)) {
			cerr << "Sfc index " << sfc_index << " wasn't cached serially" << endl;
			exit(EXIT_FAILURE);
		}

		boost::array<unsigned int, 3> indices = mapping.get_indices(sfc_index);
		serial_file << indices[0] << "\t"
			<< indices[1] << "\t"
			<< indices[2] << "\t"
			<< sfc_index
			<< endl;
	}
	serial_file.close();*/
	mapping.clear();


	// cache grid <-> sfc mapping in parallel using several Sfc instances
	vector<Sfc<3> > mappings;
	// thread i caches from index_rages[i].first to index_ranges[i].second
	vector<pair<unsigned int, unsigned int> > index_ranges;

	const time_t before_parallel = time(NULL);
	int threads = 1;

	#pragma omp parallel
	{
		#pragma omp master
		threads = omp_get_num_threads();
	}

	const unsigned int total_indices = indices[0] * indices[1] * indices[2];
	const unsigned int indices_per_thread = total_indices / threads + 1;

	for (int i = 0; i < threads; i++) {
		mappings.push_back(Sfc<3>(indices));

		// calculate which indices to cache by which thread
		unsigned int start_index = i * indices_per_thread;
		if (start_index >= total_indices) {
			start_index = total_indices - 1;
		}

		unsigned int end_index = (i + 1) * indices_per_thread;
		if (end_index >= total_indices) {
			end_index = total_indices - 1;
		}

		index_ranges.push_back(make_pair(start_index, end_index));
	}

	// cache mappings in parallel
	#pragma omp parallel for
	for (int i = 0; i < threads; i++) {
		mappings[i].cache_sfc_index_range(index_ranges[i].first, index_ranges[i].second);
	}

	const time_t after_parallel = time(NULL);

	cout << "Parallel mapping calculated in " << after_parallel - before_parallel << " seconds using " << threads << " thread(s)" << endl;

	// write out sfc traversal order calculated in parallel for verification
	/*ofstream parallel_file("parallel.dat");

	int current_mapping = 0;
	unsigned int sfc_index = index_ranges[current_mapping].first;
	while (sfc_index <= index_ranges[threads - 1].second) {

		if (!mappings[current_mapping].is_cached(sfc_index)) {
			cerr << "Sfc index " << sfc_index << " wasn't cached in parallel" << endl;
			exit(EXIT_FAILURE);
		}

		boost::array<unsigned int, 3> indices = mappings[current_mapping].get_indices(sfc_index);
		parallel_file << indices[0] << "\t"
			<< indices[1] << "\t"
			<< indices[2] << "\t"
			<< sfc_index
			<< endl;

		sfc_index++;

		if (sfc_index > index_ranges[current_mapping].second) {
			current_mapping++;
		}
	}
	parallel_file.close();*/

	return EXIT_SUCCESS;
}

