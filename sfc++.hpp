/*!
A class for mapping 2 and 3 dimensional data into 1 dimension.

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


Data is mapped so that neighboring cells in > 1d are also quite close
neighbors in 1d. Based on the following paper:
Campbell, Devine, Flaherty, Gervasio, Teresco:
Dynamic Octree Load Balancing Using Space-Filling Curves
available at http://j.teresco.org//research/publications/octpart02/octpart02.pdf

Requires http://www.boost.org/
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
A class for internal variables and functions with trivial or no dimensional dependence.
*/
template<unsigned int Dimensions> class Sfc_internal
{

public:

	/*!
	Clears stored mappings of cells in user's grid to sfc indices.
	*/
	void clear_cache(void)
	{
		assert(this->initialized);
		this->user_indices.clear();
	}


	/*!
	Returns the index in 1 dimension of a cell with given indinces.

	Indices start at 0.
	*/
	unsigned int get_index(const boost::array<unsigned int, Dimensions>& given_indices)
	{
		assert(this->initialized);

		if (this->user_indices.count(given_indices) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " ";
			for (unsigned int i = 0; i < Dimensions; i++) {
				std::cerr << given_indices[i] << " ";
			}
			std::cerr << std::endl;
			assert(false);
			//this->calculate_index(given_indices);
		}

		return this->user_indices.at(given_indices);
	}


	/*!
	Returns internal sfc coordinate traversal order.

	Might be empty.
	*/
	const std::list<std::vector<unsigned int> >& get_coordinates(void)
	{
		return this->coordinates;
	}



protected:

	bool initialized;

	// length of user's grid in indices
	boost::array<unsigned int, Dimensions> length;
	// smallest power of 2 at least as large any of above
	unsigned int uncropped_length;

	/*
	At what refinement level one sfc coordinate
	corresponds to one cell in user's grid.
	*/
	unsigned int refinement_level;

	// sfc coordinates in the order as they are traversed
	std::list<std::vector<unsigned int> > coordinates;

	// calculated mappings between cell indeices of user's grid and 1d sfc indices
	boost::unordered_map<boost::array<unsigned int, Dimensions>, unsigned int> user_indices;
	//boost::unordered_map<unsigned int, boost::unordered_map<unsigned int, unsigned int> > user_indices;


	// fille order and orientation arrays
	virtual void fill_internal_arrays(void) = 0;


	// turns one sfc coordinate into 2^d where d is the number of dimensions
	virtual void refine(std::list<std::vector<unsigned int> >::iterator position) = 0;


	// refines all sfc coordinates once
	void refine_all(void)
	{
		for (std::list<std::vector<unsigned int> >::iterator
			item = this->coordinates.begin();
			item != this->coordinates.end();
		) {
			// the last one might get removed
			std::list<std::vector<unsigned int> >::iterator to_refine = item;
			item++;
			this->refine(to_refine);
		}
	}


	// calculates and stores the mapped 1d index of all cell in the user's grid
	virtual void cache_all(void) = 0;


	// prepares the class for use
	void initialize(const boost::array<unsigned int, Dimensions>& given_length)
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
		unsigned int longest_dimension = 0;
		for (unsigned int i = 0; i < Dimensions; i++) {
			longest_dimension = std::max(longest_dimension, this->length[i]);
		}

		while ((1u << this->refinement_level) < longest_dimension) {
			this->refinement_level++;
		}

		this->uncropped_length = (1 << this->refinement_level);

		this->initialized = true;
		this->cache_all();
	}
};


/*!
The general version doesn't do anything.

Patches implementing this for N dimensions more than welcome :)
*/
template<unsigned int Dimensions> class Sfc : public Sfc_internal<Dimensions> {};


/*!
Class for mapping the cells of 2-dimensional grid to 1 dimension.

Neighboring cells in 2d are also quite close in 1d.
They are closest when the given grid is of equal length in
cells in each dimension and a power of 2.
*/
template<> class Sfc<2> : public Sfc_internal<2>
{

public:

	Sfc(const boost::array<unsigned int, 2>& given_length)
	{
		this->fill_internal_arrays();
		this->initialize(given_length);
	}


	void cache_all(void)
	{
		assert(this->initialized);
		boost::array<unsigned int, 2> indices;

		// create all sfc coordinates
		for (unsigned int i = 0; i < this->refinement_level; i++) {
			this->refine_all();
		}

		// initialize user_indices
		// TODO generalize to N dimensions
		for (unsigned int y = 0; y < this->length[1]; y++) {
			indices[0] = y;
			for (unsigned int x = 0; x < this->length[0]; x++) {
				// TODO switch x and y
				indices[1] = x;
				this->user_indices[indices] = std::numeric_limits<unsigned int>::max();
			}
		}

		// get the sfc coordinate of all cells in user's grid
		boost::unordered_map<
			std::vector<unsigned int>,
			boost::array<unsigned int, 2>
		> coord_to_indices;

		for (unsigned int y = 0; y < this->length[1]; y++) {
			indices[1] = y;
			for (unsigned int x = 0; x < this->length[0]; x++) {
				indices[0] = x;
				coord_to_indices[this->get_coordinate(indices)] = indices;
			}
		}

		// fill out the order in which sfc traverses the 2d indices
		int sfc_index = 0;
		for (std::list<std::vector<unsigned int> >::const_iterator
			coordinate = this->coordinates.begin();
			coordinate != this->coordinates.end();
			coordinate++
		) {
			if (coord_to_indices.count(*coordinate) == 0) {
				continue;
			}

			this->user_indices[coord_to_indices.at(*coordinate)] = sfc_index;
			sfc_index++;
		}
	}


	/*!
	Returns the indices in user's grid of given sfc coordinate.
	*/
	boost::array<unsigned int, 2> get_indices(const std::vector<unsigned int>& coordinate)
	{
		// current location
		boost::array<unsigned int, 2> location = boost::assign::list_of(0)(0);
		unsigned int location_diff = this->uncropped_length / 2;

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



private:

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
	Turns the sfc coordinate at given position into 4 or less longer coordinates.

	Removes all new coordinates which are outside of the user's grid.
	*/
	void refine(std::list<std::vector<unsigned int> >::iterator position)
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

		// remove those outside of the user's grid
		for (i = 0; i < 4; i++) {
			boost::array<unsigned int, 2> indices = this->get_indices(*position);
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
				position++;
			}
		}
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
	std::vector<unsigned int> get_coordinate(const boost::array<unsigned int, 2>& indices) const
	{
		std::vector<unsigned int> coordinate;

		// current search area
		boost::array<unsigned int, 2> start, end;
		for (unsigned int i = 0; i < 2; i++) {
			start[i] = 0;
			end[i] = this->uncropped_length - 1;
		}

		for (unsigned int ref_lvl = 0; ref_lvl < this->refinement_level; ref_lvl++) {
			// TODO: generalize to N dimensions
			boost::array<unsigned int, 2> before_middle;
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
Class for mapping the cells of 3-dimensional grid to 1 dimension.

Neighboring cells in 3d are also quite close in 1d.
They are closest when the given grid is of equal length in
cells in each dimension and a power of 2.
*/
template<> class Sfc<3> : public Sfc_internal<3>
{

public:

	Sfc(const boost::array<unsigned int, 3>& given_length)
	{
		this->fill_internal_arrays();
		this->initialize(given_length);
	}


	void cache_all(void)
	{
		assert(this->initialized);
		boost::array<unsigned int, 3> indices;

		// create all sfc coordinates
		for (unsigned int i = 0; i < this->refinement_level; i++) {
			this->refine_all();
		}

		// initialize user_indices
		// TODO generalize to N dimensions
		for (unsigned int z = 0; z < this->length[2]; z++) {
			indices[0] = z;
			for (unsigned int y = 0; y < this->length[1]; y++) {
				indices[1] = y;
				for (unsigned int x = 0; x < this->length[0]; x++) {
					indices[2] = x;
					this->user_indices[indices] = std::numeric_limits<unsigned int>::max();
				}
			}
		}

		// get the sfc coordinate of all cells in user's grid
		boost::unordered_map<
			std::vector<unsigned int>,
			boost::array<unsigned int, 3>
		> coord_to_indices;

		for (unsigned int z = 0; z < this->length[2]; z++) {
			indices[2] = z;
			for (unsigned int y = 0; y < this->length[1]; y++) {
				indices[1] = y;
				for (unsigned int x = 0; x < this->length[0]; x++) {
					indices[0] = x;
					coord_to_indices[this->get_coordinate(indices)] = indices;
				}
			}
		}

		// fill out the order in which sfc traverses the 2d indices
		int sfc_index = 0;
		for (std::list<std::vector<unsigned int> >::const_iterator
			coordinate = this->coordinates.begin();
			coordinate != this->coordinates.end();
			coordinate++
		) {
			if (coord_to_indices.count(*coordinate) == 0) {
				continue;
			}

			this->user_indices[coord_to_indices.at(*coordinate)] = sfc_index;
			sfc_index++;
		}
	}


	/*!
	Returns the indices in user's grid of given sfc coordinate.
	*/
	boost::array<unsigned int, 3> get_indices(const std::vector<unsigned int>& coordinate)
	{
		// current location
		boost::array<unsigned int, 3> location = boost::assign::list_of(0)(0)(0);
		unsigned int location_diff = this->uncropped_length / 2;

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



private:

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
	Turns the sfc coordinate at given position into 8 or less longer coordinates.

	Removes all new coordinates which are outside of the user's grid.
	*/
	void refine(std::list<std::vector<unsigned int> >::iterator position)
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

		// remove those outside of the user's grid
		for (i = 0; i < 8; i++) {
			boost::array<unsigned int, 3> indices = this->get_indices(*position);
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
				position++;
			}
		}
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
	std::vector<unsigned int> get_coordinate(const boost::array<unsigned int, 3>& indices) const
	{
		std::vector<unsigned int> coordinate;

		// current search area
		boost::array<unsigned int, 3> start, end;
		for (unsigned int i = 0; i < 3; i++) {
			start[i] = 0;
			end[i] = this->uncropped_length - 1;
		}

		for (unsigned int ref_lvl = 0; ref_lvl < this->refinement_level; ref_lvl++) {
			// TODO: generalize to N dimensions
			boost::array<unsigned int, 3> before_middle;
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

