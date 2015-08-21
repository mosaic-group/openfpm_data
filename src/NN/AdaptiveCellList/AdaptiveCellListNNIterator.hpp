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

// TODO: copied from CellNNIterator (including the fixed +1 at SYM!), doesn't look like a very good idea to me
#define FULL openfpm::math::pow(3,dim)
#define SYM openfpm::math::pow(3,dim)/2+1
#define CRS openfpm::math::pow(2,dim)

template<unsigned int dim, typename CellS, unsigned int NNc_size_original, unsigned int impl> class AdaptiveCellNNIterator
{

	CellS & cl;

	Point<dim+1, typename CellS::value_type> p;

	size_t currentlevel;
	
	// With FULL it always stays full, but with SYM we have to change it to FULL after the first level!
	unsigned int NNc_size;
	
	std::vector<size_t> cellindices;
	std::vector<size_t>::iterator cell_iter, cell_iter_end;
	decltype(cl.getCellContents(0).first) ele_iter, ele_iter_end;
	
	// TODO: theres got to be a more intelligent way of doing this.
	// Also: we pass NNc as a parameter, to use FULL inside SYM. See code below.
	inline void populateIndices(int current_dim, Point<dim+1, typename CellS::value_type>& p, typename CellS::value_type edgelength, unsigned int NNc) {
		if (current_dim < 0) {
			if(cl.isInside(p))
				cellindices.push_back(cl.findCellIndex(p));
			//std::cout << p.toString() << std::endl;
		}
		else if(NNc == SYM) {
			// This means that the first cell we start in is the cell of the origin particle itself!
			populateIndices(current_dim - 1, p, edgelength, SYM);
			p.get(current_dim) += edgelength;
			populateIndices(current_dim - 1, p, edgelength, FULL);
			p.get(current_dim) -= edgelength;
		}
		else { // Assumption: FULL
			populateIndices(current_dim - 1, p, edgelength, NNc);
			p.get(current_dim) -= edgelength;
			populateIndices(current_dim - 1, p, edgelength, NNc);
			p.get(current_dim) += edgelength+edgelength;
			populateIndices(current_dim - 1, p, edgelength, NNc);
			p.get(current_dim) -= edgelength;
		}
	}
	
	inline void updateCellIterators() {
		cellindices.resize(0);
		cellindices.reserve(NNc_size);
		
		Point<dim+1, typename CellS::value_type> center(p);
		center.get(dim) = cl.getRadiusForLevel(currentlevel);
		typename CellS::value_type edgelength = cl.getEdgeLengthForLevel(currentlevel);
		
		populateIndices(dim - 1, center, edgelength, NNc_size);
		
		/*
		if(NNc_size == SYM) {
			std::cout << " -- On level " << currentlevel << ": ";
			for (auto i : cellindices)
				std::cout << i << ", ";
			std::cout << std::endl;
		}
		//*/
		
		cell_iter = cellindices.begin();
		cell_iter_end = cellindices.end();
	}
	
	inline void updateEleIterators() {
		//if(NNc_size == SYM) std::cout << "Request ele_iters for cell " << *cell_iter << std::endl;
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
	AdaptiveCellNNIterator(std::pair<Point<dim+1, typename CellS::value_type>, size_t> pointpair, CellS & cl)
		:cl(cl), p(pointpair.first), currentlevel(0), NNc_size(NNc_size_original)
	{
		if(NNc_size != FULL)
			currentlevel = cl.getLevelOfRadius(p.get(dim));
		
		//std::cout << "NN of " << pointpair.second << " at " << p.toString() << " starting on lv " << currentlevel << std::endl;
		
		updateCellIterators();
		updateEleIterators();
		
		if(NNc_size == SYM) {
			// We know this cell contains an element - the original particle.
			assert(ele_iter != ele_iter_end);
			//for(auto it = ele_iter; it != ele_iter_end; ++it) std::cout << it->second << " at " << it->first.toString() << std::endl;
			// So progress to it to ensure uniqueness for SYM.
			//std::cout << "> ";
			while (ele_iter->second != pointpair.second) {
				++ele_iter;
				//std::cout << "+";
			}
			++ele_iter;
			//std::cout << " <, now at " << ele_iter->second << std::endl;
		}
		
		// Even in SYM the original particle might be the only one existing in the cell!
		while(ele_iter == ele_iter_end && currentlevel <= cl.maxlevel()) {
			++cell_iter;
			
			if(cell_iter == cell_iter_end) {
				++currentlevel;
				if(NNc_size_original == SYM)
					NNc_size = FULL;
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
		//if(NNc_size == SYM) std::cout << currentlevel << " <= " << cl.maxlevel() << std::endl;
		return currentlevel <= cl.maxlevel();
	}

	/*! \brief take the next element
	 *
	 */
	AdaptiveCellNNIterator & operator++()
	{
		// Condition: don't start in an empty cell, otherwise this `++` will trick the following `while`.
		++ele_iter;
		
		while(ele_iter == ele_iter_end && currentlevel <= cl.maxlevel()) {
			++cell_iter;
			
			if(cell_iter == cell_iter_end) {
				++currentlevel;
				if(NNc_size_original == SYM)
					NNc_size = FULL;
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
