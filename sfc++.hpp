/*
A class for mapping 2 and 3 dimensional data into 1 dimension.

Copyright 2010, 2011 Ilja Honkonen
Copyright 2011 Finnish Meteorological Institute

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

/*!
\mainpage sfc++

\section intro_sec Introduction
A class for mapping 2 and 3 dimensional data into 1 dimension.

Data is mapped so that neighboring cells in > 1d are also quite close
neighbors in 1d. Based on the following paper:
Campbell, Devine, Flaherty, Gervasio, Teresco:
Dynamic Octree Load Balancing Using Space-Filling Curves
available at http://j.teresco.org//research/publications/octpart02/octpart02.pdf

Requires http://www.boost.org/
*/

/*
TODO: See if this would speed things up:
A new algorithm for encoding and decoding the Hilbert order,
Ningtao Chen, Nengchao Wang and Baochang Shi,
Softw. Pract. Exper. 37:897–908 2007,
DOI: 10.1002/spe.793
*/

#ifndef SFC_HPP
#define SFC_HPP

#include "algorithm"
#include "boost/array.hpp"
#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "boost/functional/hash.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cstdio"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "list"
#include "map"
#include "utility"
#include "vector"

// hashing for boost arrays, from http://lists.boost.org/boost-users/2010/08/61836.php
namespace boost {
	template <typename T, std::size_t N> std::size_t hash_value(const boost::array<T, N>& array)
	{
		return boost::hash_range(array.begin(), array.end());
	}
}


namespace sfc {

/*!
A class for Sfc variables and functions with trivial or no dimensional dependence.

Doesn't implement all functions so use the Sfc class in 2d or 3d instead.
Functions common to both dimensions are documented here though.
*/
template<unsigned int Dimensions, class Index_T = unsigned int> class Sfc_common
{

public:

	/*!
	Caches all mappings between grid indices and sfc indices.

	Use this if most of the sfc is going to be traversed and memory
	isn't a limiting factor.

	Requires O(N * log2(L)) operations where N is the number of
	cells in the grid and L is the longest length of the grid
	in indices between all dimensions.
	*/
	void cache_all(void)
	{
		assert(this->initialized);

		this->cache_sfc_index_range(0, this->size() - 1);
	};


	/*!
	Caches mappings between grid indices and given sfc index range inclusive.

	Use this if only a certain range of the sfc is going to be
	traversed or cacle_all would require too much memory.

	Requires O(N * log2(L)) operations where N is the number of
	indices in given range and L is the longest length of the grid
	in indices between all dimensions.
	*/
	void cache_sfc_index_range(
		const Index_T index_start,
		const Index_T index_end
	) {
		assert(this->initialized);

		if (index_start > index_end) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "index_start must be <= index_end"
				<< std::endl;
			abort();
		}

		this->calculate_indices_range(index_start, index_end);
	};


	/*!
	Removes all stored information.

	Removes mappings between grid indices and sfc indices
	and calculated sfc coordinates.
	*/
	void clear(void)
	{
		this->user_indices.clear();
		this->sfc_indices.clear();
		this->coordinates.clear();
		this->coordinates.insert(this->coordinates.begin(), std::vector<unsigned int>());
	}


	/*!
	Removes the mapping between given indices and an sfc index from the cache.
	*/
	void clear(const boost::array<Index_T, Dimensions>& given_indices)
	{
		this->user_indices.erase(given_indices);
	}


	/*!
	Removes the mapping between given sfc index and grid indices from the cache.
	*/
	void clear(const Index_T given_index)
	{
		this->sfc_indices.erase(given_index);
	}


	/*!
	Returns the size of the grid in cells.
	*/
	Index_T size(void) const
	{
		Index_T result = 1;
		for (Index_T dimension = 0; dimension < Dimensions; dimension++) {
			result *= this->length[dimension];
		}
		return result;
	}


	/*!
	Returns the sfc index of given indices in the grid.

	Given the indices for a cell in the grid returns the corresponding index in the
	space-filling curve. For example the following figure shows the grid indices
	and the corresponding sfc indices inside the cells of a 2 by 2 grid:
\verbatim
y

^
|
 -----
1|1|0|
 -----
0|2|3|
 -----
  0 1 -> x
\endverbatim

	Requires O(2N * log2(L)) operations where N is the number of sfc coordinates
	(~cached mappings of grid indices to sfc indices) and and L is the longest
	length of the grid in indices between all dimensions.

	The returned value is also cached for faster access later.
	Indices start at 0.
	*/
	Index_T get_sfc_index(const boost::array<Index_T, Dimensions>& given_indices)
	{
		assert(this->initialized);

		if (this->user_indices.count(given_indices) == 0) {
			this->calculate_sfc_index(given_indices);
		}

		return this->user_indices.at(given_indices);
	}


	/*!
	Returns the indices in the grid of given sfc index.

	The returned value is also cached for faster access later.
	Indices start at 0.
	*/
	boost::array<Index_T, Dimensions> get_indices(const Index_T given_index)
	{
		assert(this->initialized);

		if (this->sfc_indices.count(given_index) == 0) {
			this->calculate_indices(given_index);
		}

		return this->sfc_indices.at(given_index);
	}


	/*!
	Adds cached sfc and grid indices from rhs to lhs.

	Does nothing if in both instances of Sfc_common the grid isn't identical.
	*/
	Sfc_common& operator += (const Sfc_common& rhs)
	{
		for (unsigned int i = 0; i < Dimensions; i++) {
			if (this->length[i] != rhs.length[i]) {
				return *this;
			}
		}

		// add cached grid indices
		for (typename boost::unordered_map<boost::array<Index_T, Dimensions>, Index_T>::const_iterator
			item = rhs.user_indices.begin();
			item != rhs.user_indices.end();
			item++
		) {
			if (this->user_indices.count(item->first) > 0) {
				continue;
			}

			this->user_indices[item->first] = item->second;
		}

		// add cached sfc indices
		for (typename boost::unordered_map<Index_T, boost::array<Index_T, Dimensions> >::const_iterator
			item = rhs.sfc_indices.begin();
			item != rhs.sfc_indices.end();
			item++
		) {
			if (this->sfc_indices.count(item->first) > 0) {
				continue;
			}

			this->sfc_indices[item->first] = item->second;
		}

		return *this;
	}


	/*!
	Returns true if given sfc index has been cached, false otherwise.
	*/
	bool is_cached(const Index_T given_index)
	{
		if (this->sfc_indices.count(given_index) > 0) {
			return true;
		} else {
			return false;
		}
	}



protected:

	bool initialized;

	// length of user's grid in indices
	boost::array<Index_T, Dimensions> length;
	// smallest power of 2 at least as large any of above
	Index_T uncropped_length;

	/*
	At what refinement level one sfc coordinate
	corresponds to one cell in user's grid.
	*/
	Index_T refinement_level;

	// sfc coordinates in the order as they are traversed
	std::list<std::vector<unsigned int> > coordinates;

	// calculated mappings between cell indices of user's grid and sfc indices
	boost::unordered_map<boost::array<Index_T, Dimensions>, Index_T> user_indices;

	// calculated mappings between sfc indices and cell indices of the user's grid
	boost::unordered_map<Index_T, boost::array<Index_T, Dimensions> > sfc_indices;
	// TODO use boost::bimap for the above two?


	// fille order and orientation arrays
	virtual void fill_internal_arrays(void) = 0;


	// turns one sfc coordinate into 2^d where d is the number of dimensions
	virtual std::list<std::vector<unsigned int> >::iterator refine(std::list<std::vector<unsigned int> >::iterator position) = 0;


	/*!
	Returns the indices in user's grid of given sfc coordinate.

	If the given coordinate spans more than one set of indices
	the one closest to (0, 0) is returned.
	*/
	virtual boost::array<Index_T, Dimensions>
	get_indices_from_coordinate(const std::vector<unsigned int>& coordinate) = 0;


	/*!
	Sets the unropped length of the grid, sets the initial coordinate, etc.
	*/
	void initialize(const boost::array<Index_T, Dimensions>& given_length)
	{
		for (unsigned int i = 0; i < Dimensions; i++) {
			if (given_length[i] == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Zero length in dimension " << i
					<< std::endl;
				abort();
			}

			this->length[i] = given_length[i];
		}

		// the first coordinate is empty
		this->coordinates.insert(this->coordinates.begin(), std::vector<unsigned int>());

		// get first power of two not smaller than any given dimension
		this->refinement_level = 0;
		Index_T longest_dimension = 0;
		for (unsigned int i = 0; i < Dimensions; i++) {
			longest_dimension = std::max(longest_dimension, this->length[i]);
		}

		while ((Index_T(1) << this->refinement_level) < longest_dimension) {
			this->refinement_level++;
		}

		this->uncropped_length = (1 << this->refinement_level);

		this->initialized = true;
	}


	// refines all sfc coordinates once
	void refine_all(void)
	{
		for (std::list<std::vector<unsigned int> >::iterator
			item = this->coordinates.begin();
			item != this->coordinates.end();
		) {
			item = this->refine(item);
			item++;
		}
	}


	/*!
	Returns the 1d length of given coordinate in indices of the user's grid.
	*/
	Index_T get_length_in_user_indices(const std::vector<unsigned int>& coordinate) const
	{
		return this->uncropped_length / (Index_T(1) << coordinate.size());
	}


	/*!
	Returns the number of cells in the user's grid that the given coordinate encompases.
	*/
	Index_T get_number_of_user_indices(const std::vector<unsigned int>& coordinate)
	{
		const Index_T length = this->get_length_in_user_indices(coordinate);
		if (length == 0 || length > this->uncropped_length) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Invalid length in indices: " << length
				<< std::endl;
			abort();
		}

		const boost::array<Index_T, Dimensions>& indices = this->get_indices_from_coordinate(coordinate);

		Index_T number_of = 1;
		for (unsigned int dimension = 0; dimension < Dimensions; dimension++) {
			if (this->length[dimension] <= indices[dimension]) {
				return 0;
			}

			number_of *= std::min(length, this->length[dimension] - indices[dimension]);
		}

		return number_of;
	}


	/*!
	Adds indices in the grid corresponding to given sfc index into the cache.

	Returns the corresponding sfc coordinate.
	*/
	std::list<std::vector<unsigned int> >::iterator calculate_indices(const Index_T given_index)
	{
		if (given_index >= this->size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Given index is outside of the grid: " << given_index
				<< std::endl;
			abort();
		}

		std::list<std::vector<unsigned int> >::iterator position = this->coordinates.begin();
		Index_T length = this->get_number_of_user_indices(*position);
		Index_T sfc_index = 0;

		// find the coordinate with given sfc index inside
		while (!(sfc_index <= given_index
		&& sfc_index + length > given_index)) {
			sfc_index += length;
			// FIXME: goes out of the list if given_index == last index
			position++;
			length = this->get_number_of_user_indices(*position);
		}

		// refine the coordinate until it represents only one sfc index
		while (length > 1) {
			sfc_index += length;
			position = this->refine(position);
			length = this->get_number_of_user_indices(*position);
			sfc_index -= length;

			// find the coordinate with given indices inside			
			while (!(sfc_index <= given_index
			&& sfc_index + length > given_index)) {
				position--;
				length = this->get_number_of_user_indices(*position);
				sfc_index -= length;
			}
		}

		// calculate indices in the user's grid corresponding to given sfc index
		const boost::array<Index_T, Dimensions>& indices = this->get_indices_from_coordinate(*position);
		this->sfc_indices[sfc_index] = indices;
		this->user_indices[indices] = sfc_index;

		return position;
	}


	void calculate_indices_range(
		const Index_T index_start,
		const Index_T index_end
	) {
		if (index_start > index_end) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "index_start must be <= index_end"
				<< std::endl;
			abort();
		}

		const std::list<std::vector<unsigned int> >::iterator start_position = this->calculate_indices(index_start);

		if (index_start == index_end) {
			return;
		}

		// calculate last requested sfc index and all in between
		Index_T sfc_index = index_end;
		std::list<std::vector<unsigned int> >::iterator position = this->calculate_indices(index_end);

		position--;
		while (position != start_position) {
			Index_T length = this->get_number_of_user_indices(*position);
			while (length > 1) {
				position = this->refine(position);
				length = this->get_number_of_user_indices(*position);
			}
			sfc_index--;

			const boost::array<Index_T, Dimensions>& indices = this->get_indices_from_coordinate(*position);

			this->sfc_indices[sfc_index] = indices;
			this->user_indices[indices] = sfc_index;

			position--;
		}
	}


	/*!
	Adds the sfc index of given indices in the user's grid into the cache.
	*/
	void calculate_sfc_index(const boost::array<Index_T, Dimensions>& given_indices)
	{
		for (unsigned int dimension = 0; dimension < Dimensions; dimension++) {
			if (given_indices[dimension] >= this->length[dimension]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Invalid index in dimension " << dimension
					<< ": " << given_indices[dimension]
					<< std::endl;
				abort();
			}
		}

		std::list<std::vector<unsigned int> >::iterator position = this->coordinates.begin();
		Index_T length = this->get_length_in_user_indices(*position);
		boost::array<Index_T, Dimensions> indices = this->get_indices_from_coordinate(*position);

		// find the coordinate with given indices inside
		while (!this->cells_overlap(indices, length, given_indices)) {
			position++;
			length = this->get_length_in_user_indices(*position);
			indices = this->get_indices_from_coordinate(*position);
		}

		// refine the coordinate until it represents only one set of user indices
		while (length > 1) {
			position = this->refine(position);
			length = this->get_length_in_user_indices(*position);
			indices = this->get_indices_from_coordinate(*position);

			// find the coordinate with given indices inside			
			while (!this->cells_overlap(indices, length, given_indices)) {
				position--;
				indices = this->get_indices_from_coordinate(*position);
			}
		}

		// calculate the sfc index of found coordinate
		int sfc_index = 0;
		std::list<std::vector<unsigned int> >::iterator temp = this->coordinates.begin();
		while (temp != position) {
			sfc_index += get_number_of_user_indices(*temp);
			temp++;
		}

		this->user_indices[given_indices] = sfc_index;
		this->sfc_indices[sfc_index] = given_indices;
	}


	/*!
	Returns true if cell2 is inside cell1.

	indices1 is the starting point of cell1 in the grid (closest to origin).
	length1 is the length of cell1 in indices in every dimension.
	The length in indices of cell2 is assumed to be 1.
	*/
	bool cells_overlap(
		const boost::array<Index_T, Dimensions>& indices1,
		const Index_T length1,
		const boost::array<Index_T, Dimensions>& indices2
	) const
	{
		for (unsigned int dimension = 0; dimension < Dimensions; dimension++) {
			if (!(indices1[dimension] <= indices2[dimension]
			&& indices1[dimension] + length1 > indices2[dimension])) {
				return false;
			}
		}

		return true;
	}
};


/*!
A Sfc class for an arbitrary number of dimensions, doesn't do anything.

See documentation of Sfc<2> and Sfc<3> specializations.
Patches implementing this for N dimensions more than welcome :)
*/
template<unsigned int Dimensions, class Index_T = unsigned int> class Sfc : public Sfc_common<Dimensions, Index_T> {};


/*!
Class for mapping the cells of a 2-dimensional grid to 1 dimension.
*/
template<class Index_T> class Sfc<2, Index_T> : public Sfc_common<2, Index_T>
{

public:

	/*!
	Initializes a sfc mapping for a 2d grid of given size in indices.

	Best results (e.g. neighboring cells in 2d are closest in 1d) are
	obtained if given lengths are powers of 2. Results are probably
	also as good if given lengths are multiples of 2. Otherwise some
	cells' neighbors in 2d might be a few cells further away than
	usual in 1d.

	If the whole sfc is going to be traversed and memory is not a
	limiting factor the fastest way to calculate all the mappings
	is to use the cache_all method.
	*/
	Sfc(const boost::array<Index_T, 2>& given_length)
	{
		// TODO possible to move this into Sfc_common?
		this->fill_internal_arrays();
		this->initialize(given_length);
	}



private:

	/*
	The logic for using these is described in
	http://j.teresco.org//research/publications/octpart02/octpart02.pdf
	*/
	boost::array<boost::array<unsigned int, 4>, 4> order;
	boost::array<boost::array<unsigned int, 4>, 4> orientation;

	void fill_internal_arrays(void)
	{
		this->order = boost::assign::list_of
			(boost::assign::list_of(0u)(1u)(3u)(2u))
			(boost::assign::list_of(0u)(2u)(3u)(1u))
			(boost::assign::list_of(3u)(1u)(0u)(2u))
			(boost::assign::list_of(3u)(2u)(0u)(1u));

		this->orientation = boost::assign::list_of
			(boost::assign::list_of(1u)(0u)(0u)(2u))
			(boost::assign::list_of(0u)(1u)(1u)(3u))
			(boost::assign::list_of(3u)(2u)(2u)(0u))
			(boost::assign::list_of(2u)(3u)(3u)(1u));
	}


	/*!
	Returns the indices in the grid of given sfc coordinate.
	*/
	boost::array<Index_T, 2> get_indices_from_coordinate(const std::vector<unsigned int>& coordinate)
	{
		// current location
		boost::array<Index_T, 2> location = boost::assign::list_of(0)(0);
		Index_T location_diff = this->uncropped_length / 2;

		BOOST_FOREACH(unsigned int number, coordinate) {
			switch (number) {
			case 0:
				location[0] += location_diff;
				location[1] += location_diff;
				break;
			case 1:
				location[1] += location_diff;
				break;
			case 2:
				location[0] += location_diff;
				break;
			case 3:
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Invalid number in coordinate: " << number
					<< std::endl;
				abort();
				break;
			}

			location_diff /= 2;
		}

		return location;
	}


	/*!
	Turns the sfc coordinate at given position into 4 or less longer coordinates.

	Removes all new coordinates which are outside of the user's grid.
	Returns an iterator to the last created coordinate.
	*/
	std::list<std::vector<unsigned int> >::iterator refine(std::list<std::vector<unsigned int> >::iterator position)
	{
		int orientation = this->get_orientation(*position);

		// append offsprings in correct order
		this->coordinates.insert(position, 3, *position);
		unsigned int i = 3;
		while (true) {
			position->push_back(this->order[orientation][i]);
			position--;
			if (i == 0) {
				break;
			} else {
				i--;
			}

		}
		position++;

		std::list<std::vector<unsigned int> >::iterator last_created = position;

		// remove those outside of the user's grid
		for (i = 0; i < 4; i++) {
			boost::array<Index_T, 2> indices = this->get_indices_from_coordinate(*position);
			bool is_outside = false;
			for (unsigned int dimension = 0; dimension < 2; dimension++) {
				if (indices[dimension] >= this->length[dimension]) {
					is_outside = true;
					break;
				}
			}

			if (is_outside) {
				std::list<std::vector<unsigned int> >::iterator to_erase = position;
				position++;
				this->coordinates.erase(to_erase);
			} else {
				last_created = position;
				position++;
			}
		}

		return last_created;
	}


	/*!
	Returns the orientation of a space-filling curve coordinate.
	*/
	unsigned int get_orientation(const std::vector<unsigned int>& coordinate)
	{
		unsigned int orientation = 0;

		for (std::vector<unsigned int>::const_iterator
			number = coordinate.begin();
			number != coordinate.end();
			number++
		) {
			unsigned int i = 0;
			while (this->order[orientation][i] != *number) {
				i++;
			}
			orientation = this->orientation[orientation][i];
		}

		return orientation;
	}


	/*!
	Returns the space-filling curve coordinate of given 2d indices.
	*/
	std::vector<unsigned int> get_coordinate(const boost::array<Index_T, 2>& indices) const
	{
		std::vector<unsigned int> coordinate;

		// current search area
		boost::array<Index_T, 2> start, end;
		for (unsigned int i = 0; i < 2; i++) {
			start[i] = 0;
			end[i] = this->uncropped_length - 1;
		}

		for (Index_T ref_lvl = 0; ref_lvl < this->refinement_level; ref_lvl++) {
			// TODO: generalize to N dimensions
			boost::array<Index_T, 2> before_middle;
			for (unsigned int i = 0; i < 2; i++) {
				before_middle[i] = start[i] + (end[i] - start[i]) / 2;
			}

			if (indices[0] <= before_middle[0]) {

				if (indices[1] <= before_middle[1]) {
					coordinate.push_back(3);
					end[0] = before_middle[0];
					end[1] = before_middle[1];
				} else {
					coordinate.push_back(1);
					end[0] = before_middle[0];
					start[1] = before_middle[1] + 1;
				}

			} else {

				if (indices[1] <= before_middle[1]) {
					coordinate.push_back(2);
					start[0] = before_middle[0] + 1;
					end[1] = before_middle[1];
				} else {
					coordinate.push_back(0);
					start[0] = before_middle[0] + 1;
					start[1] = before_middle[1] + 1;
				}
			}
		}

		return coordinate;
	}

};	// class Sfc<2>


/*!
Class for mapping the cells of a 3-dimensional grid to 1 dimension.
*/
template<class Index_T> class Sfc<3, Index_T> : public Sfc_common<3, Index_T>
{

public:

	/*!
	Initializes a sfc mapping for a 3d grid of given size in indices.

	Best results (e.g. neighboring cells in 3d are closest in 1d) are
	obtained if given lengths are powers of 2. Results are probably
	also as good if given lengths are multiples of 2. Otherwise some
	cells' neighbors in 3d might be a few cells further away than
	usual in 1d.

	If the whole sfc is going to be traversed and memory is not a
	limiting factor the fastest way to calculate all the mappings
	is to use the cache_all method.
	*/
	Sfc(const boost::array<Index_T, 3>& given_length)
	{
		// TODO possible to move this into Sfc_common?
		this->fill_internal_arrays();
		this->initialize(given_length);
	}



private:

	/*
	The logic for using these is described in
	http://j.teresco.org//research/publications/octpart02/octpart02.pdf
	*/
	boost::array<boost::array<unsigned int, 8>, 24> order;
	boost::array<boost::array<unsigned int, 8>, 24> orientation;

	void fill_internal_arrays(void)
	{
		this->order = boost::assign::list_of
			(boost::assign::list_of(0u)(1u)(3u)(2u)(6u)(7u)(5u)(4))
			(boost::assign::list_of(0u)(4u)(6u)(2u)(3u)(7u)(5u)(1))
			(boost::assign::list_of(0u)(1u)(5u)(4u)(6u)(7u)(3u)(2))
			(boost::assign::list_of(5u)(1u)(0u)(4u)(6u)(2u)(3u)(7))
			(boost::assign::list_of(3u)(7u)(6u)(2u)(0u)(4u)(5u)(1))
			(boost::assign::list_of(6u)(7u)(3u)(2u)(0u)(1u)(5u)(4))
			(boost::assign::list_of(5u)(1u)(3u)(7u)(6u)(2u)(0u)(4))
			(boost::assign::list_of(0u)(4u)(5u)(1u)(3u)(7u)(6u)(2))
			(boost::assign::list_of(5u)(4u)(0u)(1u)(3u)(2u)(6u)(7))
			(boost::assign::list_of(5u)(4u)(6u)(7u)(3u)(2u)(0u)(1))
			(boost::assign::list_of(0u)(2u)(3u)(1u)(5u)(7u)(6u)(4))
			(boost::assign::list_of(6u)(4u)(0u)(2u)(3u)(1u)(5u)(7))
			(boost::assign::list_of(5u)(7u)(3u)(1u)(0u)(2u)(6u)(4))
			(boost::assign::list_of(3u)(7u)(5u)(1u)(0u)(4u)(6u)(2))
			(boost::assign::list_of(6u)(4u)(5u)(7u)(3u)(1u)(0u)(2))
			(boost::assign::list_of(0u)(2u)(6u)(4u)(5u)(7u)(3u)(1))
			(boost::assign::list_of(6u)(2u)(0u)(4u)(5u)(1u)(3u)(7))
			(boost::assign::list_of(6u)(2u)(3u)(7u)(5u)(1u)(0u)(4))
			(boost::assign::list_of(3u)(2u)(0u)(1u)(5u)(4u)(6u)(7))
			(boost::assign::list_of(6u)(7u)(5u)(4u)(0u)(1u)(3u)(2))
			(boost::assign::list_of(5u)(7u)(6u)(4u)(0u)(2u)(3u)(1))
			(boost::assign::list_of(3u)(2u)(6u)(7u)(5u)(4u)(0u)(1))
			(boost::assign::list_of(3u)(1u)(0u)(2u)(6u)(4u)(5u)(7))
			(boost::assign::list_of(3u)(1u)(5u)(7u)(6u)(4u)(0u)(2));

		this->orientation = boost::assign::list_of
			(boost::assign::list_of( 1u)( 2u)( 0u)( 3u)( 4u)( 0u)( 5u)( 6u))
			(boost::assign::list_of( 0u)( 7u)( 1u)( 8u)( 5u)( 1u)( 4u)( 9u))
			(boost::assign::list_of(15u)( 0u)( 2u)(22u)(20u)( 2u)(19u)(23u))
			(boost::assign::list_of(20u)( 6u)( 3u)(23u)(15u)( 3u)(16u)(22u))
			(boost::assign::list_of(22u)(13u)( 4u)(12u)(11u)( 4u)( 1u)(20u))
			(boost::assign::list_of(11u)(19u)( 5u)(20u)(22u)( 5u)( 0u)(12u))
			(boost::assign::list_of( 9u)( 3u)( 6u)( 2u)(21u)( 6u)(17u)( 0u))
			(boost::assign::list_of(10u)( 1u)( 7u)(11u)(12u)( 7u)(13u)(14u))
			(boost::assign::list_of(12u)( 9u)( 8u)(14u)(10u)( 8u)(18u)(11u))
			(boost::assign::list_of( 6u)( 8u)( 9u)( 7u)(17u)( 9u)(21u)( 1u))
			(boost::assign::list_of( 7u)(15u)(10u)(16u)(13u)(10u)(12u)(17u))
			(boost::assign::list_of( 5u)(14u)(11u)( 9u)( 0u)(11u)(22u)( 8u))
			(boost::assign::list_of( 8u)(20u)(12u)(19u)(18u)(12u)(10u)( 5u))
			(boost::assign::list_of(18u)( 4u)(13u)( 5u)( 8u)(13u)( 7u)(19u))
			(boost::assign::list_of(17u)(11u)(14u)( 1u)( 6u)(14u)(23u)( 7u))
			(boost::assign::list_of( 2u)(10u)(15u)(18u)(19u)(15u)(20u)(21u))
			(boost::assign::list_of(19u)(17u)(16u)(21u)( 2u)(16u)( 3u)(18u))
			(boost::assign::list_of(14u)(16u)(17u)(15u)(23u)(17u)( 6u)(10u))
			(boost::assign::list_of(13u)(21u)(18u)(17u)( 7u)(18u)( 8u)(16u))
			(boost::assign::list_of(16u)( 5u)(19u)( 4u)( 3u)(19u)( 2u)(13u))
			(boost::assign::list_of( 3u)(12u)(20u)(13u)(16u)(20u)(15u)( 4u))
			(boost::assign::list_of(23u)(18u)(21u)(10u)(14u)(21u)( 9u)(15u))
			(boost::assign::list_of( 4u)(23u)(22u)( 6u)( 1u)(22u)(11u)( 3u))
			(boost::assign::list_of(21u)(22u)(23u)( 0u)( 9u)(23u)(14u)( 2u));
	}


	/*!
	Returns the indices in the grid of given sfc coordinate.
	*/
	boost::array<Index_T, 3> get_indices_from_coordinate(const std::vector<unsigned int>& coordinate)
	{
		// current location
		boost::array<Index_T, 3> location = boost::assign::list_of(0)(0)(0);
		Index_T location_diff = this->uncropped_length / 2;

		BOOST_FOREACH(unsigned int number, coordinate) {
			switch (number) {
			case 0:
				location[0] += location_diff;
				location[1] += location_diff;
				location[2] += location_diff;
				break;
			case 1:
				location[1] += location_diff;
				location[2] += location_diff;
				break;
			case 2:
				location[0] += location_diff;
				location[2] += location_diff;
				break;
			case 3:
				location[2] += location_diff;
				break;
			case 4:
				location[0] += location_diff;
				location[1] += location_diff;
				break;
			case 5:
				location[1] += location_diff;
				break;
			case 6:
				location[0] += location_diff;
				break;
			case 7:
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Invalid number in coordinate: " << number
					<< std::endl;
				abort();
				break;
			}

			location_diff /= 2;
		}

		return location;
	}


	/*!
	Turns the sfc coordinate at given position into 8 or less longer coordinates.

	Removes all new coordinates which are outside of the user's grid.
	*/
	std::list<std::vector<unsigned int> >::iterator refine(std::list<std::vector<unsigned int> >::iterator position)
	{
		int orientation = this->get_orientation(*position);

		// append offsprings in correct order
		this->coordinates.insert(position, 7, *position);
		unsigned int i = 7;
		while (true) {
			position->push_back(this->order[orientation][i]);
			position--;
			if (i == 0) {
				break;
			} else {
				i--;
			}

		}
		position++;

		std::list<std::vector<unsigned int> >::iterator last_created = position;

		// remove those outside of the user's grid
		for (i = 0; i < 8; i++) {
			boost::array<Index_T, 3> indices = this->get_indices_from_coordinate(*position);
			bool is_outside = false;
			for (unsigned int dimension = 0; dimension < 3; dimension++) {
				if (indices[dimension] >= this->length[dimension]) {
					is_outside = true;
					break;
				}
			}

			if (is_outside) {
				std::list<std::vector<unsigned int> >::iterator to_erase = position;
				position++;
				this->coordinates.erase(to_erase);
			} else {
				last_created = position;
				position++;
			}
		}

		return last_created;
	}


	/*!
	Returns the orientation of a space-filling curve coordinate.
	*/
	unsigned int get_orientation(const std::vector<unsigned int>& coordinate)
	{
		unsigned int orientation = 0;

		// TODO move into general version
		for (std::vector<unsigned int>::const_iterator
			number = coordinate.begin();
			number != coordinate.end();
			number++
		) {
			unsigned int i = 0;
			while (this->order[orientation][i] != *number) {
				i++;
			}
			orientation = this->orientation[orientation][i];
		}

		return orientation;
	}


	/*!
	Returns the space-filling curve coordinate of given 2d indices.
	*/
	std::vector<unsigned int> get_coordinate(const boost::array<Index_T, 3>& indices) const
	{
		std::vector<unsigned int> coordinate;

		// current search area
		boost::array<Index_T, 3> start, end;
		for (unsigned int i = 0; i < 3; i++) {
			start[i] = 0;
			end[i] = this->uncropped_length - 1;
		}

		for (Index_T ref_lvl = 0; ref_lvl < this->refinement_level; ref_lvl++) {
			// TODO: generalize to N dimensions
			boost::array<Index_T, 3> before_middle;
			for (unsigned int i = 0; i < 3; i++) {
				before_middle[i] = start[i] + (end[i] - start[i]) / 2;
			}

			if (indices[0] <= before_middle[0]) {

				if (indices[1] <= before_middle[1]) {

					if (indices[2] <= before_middle[2]) {
						coordinate.push_back(7);
						end[0] = before_middle[0];
						end[1] = before_middle[1];
						end[2] = before_middle[2];
					} else {
						coordinate.push_back(3);
						end[0] = before_middle[0];
						end[1] = before_middle[1];
						start[2] = before_middle[2] + 1;
					}

				} else {

					if (indices[2] <= before_middle[2]) {
						coordinate.push_back(5);
						end[0] = before_middle[0];
						start[1] = before_middle[1] + 1;
						end[2] = before_middle[2];
					} else {
						coordinate.push_back(1);
						end[0] = before_middle[0];
						start[1] = before_middle[1] + 1;
						start[2] = before_middle[2] + 1;
					}
				}

			} else {

				if (indices[1] <= before_middle[1]) {

					if (indices[2] <= before_middle[2]) {
						coordinate.push_back(6);
						start[0] = before_middle[0] + 1;
						end[1] = before_middle[1];
						end[2] = before_middle[2];
					} else {
						coordinate.push_back(2);
						start[0] = before_middle[0] + 1;
						end[1] = before_middle[1];
						start[2] = before_middle[2] + 1;
					}

				} else {

					if (indices[2] <= before_middle[2]) {
						coordinate.push_back(4);
						start[0] = before_middle[0] + 1;
						start[1] = before_middle[1] + 1;
						end[2] = before_middle[2];
					} else {
						coordinate.push_back(0);
						start[0] = before_middle[0] + 1;
						start[1] = before_middle[1] + 1;
						start[2] = before_middle[2] + 1;
					}
				}
			}
		}

		return coordinate;
	}

};	// class Sfc<3>

}	// namespace sfc

#endif

