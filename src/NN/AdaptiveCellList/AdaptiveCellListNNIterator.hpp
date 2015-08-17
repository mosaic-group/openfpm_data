/*
 * AdaptiveCellListNNIterator.hpp
 *
 *  Created on: May 4, 2015
 *      Author: ...
 */

#ifndef ADAPTIVECELLLISTNNITERATOR_HPP_
#define ADAPTIVECELLLISTNNITERATOR_HPP_

#include "util/mathutil.hpp"
#include <vector>
#include <iostream>

template<unsigned int dim, typename CellS, unsigned int NNc_size, unsigned int impl> class AdaptiveCellNNIterator
{

	CellS & cl;

	Point<dim+1, typename CellS::value_type> p;

	size_t currentlevel;
	
	std::vector<size_t> cellindices;
	std::vector<size_t>::iterator cell_iter, cell_iter_end;
	decltype(cl.getCellContents(0).first) ele_iter, ele_iter_end;
	
	// TODO: theres got to be a more intelligent way of doing this.
	inline void populateIndices(int current_dim, Point<dim+1, typename CellS::value_type>& p, typename CellS::value_type edgelength) {
		if (current_dim == dim) {
			if(cl.isInside(p))
				cellindices.push_back(cl.findCellIndex(p));
			//std::cout << p.toString() << std::endl;
		}
		else {
			populateIndices(current_dim + 1, p, edgelength);
			p.get(current_dim) -= edgelength;
			populateIndices(current_dim + 1, p, edgelength);
			p.get(current_dim) += edgelength+edgelength;
			populateIndices(current_dim + 1, p, edgelength);
			p.get(current_dim) -= edgelength;
		}
	}
	
	inline void updateCellIterators() {
		cellindices.resize(0);
		cellindices.reserve(openfpm::math::pow(3,dim));
		
		Point<dim+1, typename CellS::value_type> center(p);
		center.get(dim) = cl.getRadiusForLevel(currentlevel);
		typename CellS::value_type edgelength = cl.getEdgeLengthForLevel(currentlevel);
		
		populateIndices(0,center,edgelength);
		
		/*
		std::cout << " -- On level " << currentlevel << ": ";
		for (auto i : cellindices)
			std::cout << i << ", ";
		std::cout << std::endl;
		//*/
		
		cell_iter = cellindices.begin();
		cell_iter_end = cellindices.end();
	}
	
	inline void updateEleIterators() {
		//std::cout << "Request ele_iters for cell " << *cell_iter << std::endl;
		auto iterpair = cl.getCellContents(*cell_iter);
		ele_iter = iterpair.first;
		ele_iter_end = iterpair.second;
		//std::cout << "It contains " << std::distance(ele_iter, ele_iter_end) << " elements" << std::endl;
	}
	
public:
	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 */
	AdaptiveCellNNIterator(Point<dim+1, typename CellS::value_type> point, CellS & cl)
		:cl(cl), p(point), currentlevel(0)
	{
		//std::cout << "NN of " << point.toString() << std::endl;
		
		updateCellIterators();
		updateEleIterators();
		
		while(ele_iter == ele_iter_end && currentlevel <= cl.maxlevel() + 1) {
			++cell_iter;
			if(cell_iter == cell_iter_end) {
				++currentlevel;
				updateCellIterators();
			}
			updateEleIterators();
		}
		
		//std::cout << "done constructing with cell_iter at " << *cell_iter << std::endl;
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	bool isNext()
	{
		//std::cout << currentlevel << " < " << cl.maxlevel() << std::endl;
		return currentlevel <= cl.maxlevel();
	}

	/*! \brief take the next element
	 *
	 */
	AdaptiveCellNNIterator & operator++()
	{
		// Condition: don't start in an empty cell, otherwise this `++` will trick the following `if`.
		++ele_iter;
		
		while(ele_iter == ele_iter_end && currentlevel <= cl.maxlevel() + 1) {
			++cell_iter;
			if(cell_iter == cell_iter_end) {
				++currentlevel;
				updateCellIterators();
			}
			updateEleIterators();
		}
		
		return *this;
	}

	/*! \brief Get actual element
	 *
	 * \return  the actual element
	 *
	 */
	inline std::pair<Point<dim+1,typename CellS::value_type>, size_t>& get()
	{
		//std::cout << "serving " << ele_iter->first.toString() << std::endl << std::flush;
		return *ele_iter;
	}
};


#endif /* ADAPTIVECELLLISTNNITERATOR_HPP_ */
