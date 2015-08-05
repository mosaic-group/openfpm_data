#ifndef ADAPTIVECELLLIST_HPP_
#define ADAPTIVECELLLIST_HPP_

#include "../CellList/CellList.hpp"
#include "../CellList/CellDecomposer.hpp"
#include "AdaptiveCellListNNIterator.hpp"
#include "Space/SpaceBox.hpp"

#include <string>

/*
 * How does it work (normal cell lists)?
 * 
 * I have a particle and want to know all its interactions with particles in its subdomain.
 * Since I don't want to check against all other particles, I only check the proximity condition for those that are "nearby enough".
 * We could do so by dividing the subdomain into many cells of equal size (maximal cutoff/radius of the particles as the length of each side).
 * Once we accomplished that we only have to check for possible interactions with all the particles that are in neighbor cells of the cell our initially chosen particle was in.
 * This is what the NNIterator does: given some start cell list all elements elements in this start cell could maybe interact with.
 * There are multiple ways to define neighbor cells, I can check all neighbor cells (full) or only those that are "further ahaed" (sym).
 * If we know that interactions are symmetric, this makes them unique (and saves us nearly half the computations!)
 * 
 * First, let us screw all that fancy logic and implement it as naively as possible... so let us compute with all points! \o/
 */


// Stub implementation
template<unsigned int dim, typename T,  unsigned int impl=BALANCED, typename ele_container=openfpm::vector<std::pair<Point<dim+1, T>, size_t>> >
class AdaptiveCellList
{
};

/*! \brief Class for Adaptive cell list implementation
 *
 * This class implement the Adaptive cell list, ...
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam base basic object
 *
 */
template<unsigned int dim, typename T, typename ele_container>
class AdaptiveCellList<dim,T,BALANCED,ele_container> // : public CellDecomposer_sm<dim,T> // the decomposition on the toplevel stays simple for now
{
	// vector of all elements, will be sorted upon construction
	ele_container all_eles;

	//A point
	Point<dim,T> orig;
	
	// the box
	SpaceBox<dim,T> sbox;
	// minimum edge length of the box
	T D_m;
	// number of the last level
	unsigned int max_level;
	
	// AR tree represented as a hashed map from node indices to node contents (which themselves are encoded as iterators at first and last element)
	std::unordered_map<size_t, std::pair<decltype(all_eles.begin()), decltype(all_eles.end())>> cells;
	
	/*! \brief get the index (label) of the tree node at the specified position of the tree
	 * 
	 */
	inline size_t getIndex(unsigned int current_level, size_t parent_row_no, size_t current_sibling_no) {
		size_t cellindex = (parent_row_no << dim) + current_sibling_no;
		for(size_t i=0; i<current_level; i++)
			cellindex += (1l << (static_cast<size_t>(dim) * i));
		return cellindex;
	}
	
	/*! \brief Spatial sorting
	 * 
	 * \note called from construct, no need for any further calls
	 * 
	 */
	void sort(
			decltype(all_eles.begin()) begin_iter,
			decltype(all_eles.end()) end_iter,
			int current_level,
			int target_level,
			unsigned int partdim,
			SpaceBox<dim,T> box,
			size_t parent_row_no,
			size_t current_sibling_no)
	{
		//std::cout << "On lv. " << current_level << "/" << target_level << ", dim. " << partdim << ", box: ("
		//		<< box.getLow(0) << "," << box.getLow(1) << ")-(" << box.getHigh(0) << "," << box.getHigh(1) << "): "
		//		<< std::distance(begin_iter, end_iter) << std::endl;
		
		if(begin_iter == end_iter)
			return;
		
		if (partdim >= dim) { // no dimensions left to partition in (even the current one is too high)
			if(current_level < target_level) // we successfully partitioned the particles in all dimensions on the current level, let's go on and do the same thing for the next level
				sort(begin_iter, end_iter, current_level+1, target_level, 0, box, (parent_row_no << dim) + current_sibling_no, 0);
			else { //we've actually reached the target level, these are the cells we wanted to describe!
				size_t cellindex = getIndex(current_level, parent_row_no, current_sibling_no);
				// we know this part won't change, since we start sorting the vector from the beginning, so we can store iterators in our big vector
				//std::cout << "It's a " << cellindex << std::endl;
				cells.insert(std::make_pair(cellindex, std::make_pair(begin_iter, end_iter)));
			}
		}
		else { // do partition in the given dimension
			const T pivot = (box.getHigh(partdim) + box.getLow(partdim)) / 2;
			auto middle_iter = std::partition(begin_iter, end_iter,
					[&partdim, &pivot](const std::pair<Point<dim+1, T>, size_t> &a){return a.first.get(partdim) < pivot;});
			//std::cout << "With (dim. " << partdim << ") pivot " << pivot << ": "
			//		<< std::distance(begin_iter, middle_iter) << " + " << std::distance(middle_iter, end_iter) << std::endl;
			
			SpaceBox<dim,T> leftbox(box), rightbox(box);
			leftbox.setHigh(partdim, pivot);
			rightbox.setLow(partdim, pivot);
			sort(begin_iter, middle_iter, current_level, target_level, partdim+1, leftbox, parent_row_no, current_sibling_no);
			sort(middle_iter, end_iter, current_level, target_level, partdim+1, rightbox, parent_row_no, current_sibling_no + (1l << (dim-1 - partdim)));
		}
	}
	
	size_t firstIndexOfLevel(size_t level) {
		//sum i from 0 to (target_level-1): 2^(dim*i) = (2^(dim+dim*(target_level-1)) - 1) / (2^dim - 1)
		return ((1l << (dim*level)) - 1) / ((1l << dim) - 1);
	}
	
public:

	// Object type that the structure store
	typedef T value_type;

	inline size_t size() {return all_eles.size();}

	/*! Initialize the cell list
	 *
	 * \param sbox Domain where this cell list is living
	 * \param origin of the Cell list
	 *
	 */
	void Initialize(SpaceBox<dim,T> & sbox, Point<dim,T> & orig)
	{
		this->sbox = sbox;
		this->orig = orig;
		
		D_m = std::numeric_limits<T>::max();
		T tmpval;
		for(unsigned int i=0; i<dim; ++i) {
			tmpval = sbox.getHigh(i) - sbox.getLow(i);
			if(tmpval < D_m) D_m = tmpval;
		}
	}

	/*! \brief Cell list
	 *
	 * \param ... (Non default constructor if needed)
	 *
	 */
	AdaptiveCellList(SpaceBox<dim,T> & sbox, Point<dim,T> & orig)
	{
		Initialize(sbox,orig);
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate +1 (the last is the radius of the particle)
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim+1], size_t ele)
	{
		all_eles.add(std::make_pair(pos, ele));
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim+1,T> & pos, size_t ele)
	{
		all_eles.add(std::make_pair(pos, ele));
	}
	
	/*! \brief Sorts elements according so they are coherent in cells
	 * 
	 *  TODO: use some sort of isvalid-bool that is unset in add(), set in construct() and checked in the iterator
	 * 
	 * \note Has to be called after insertion and before usage!
	 * 
	 */
	void construct()
	{
		// Sort particles in descending order of their radius
		std::sort(all_eles.begin(), all_eles.end(),
				  [](std::pair<Point<dim+1, T>, size_t> a, std::pair<Point<dim+1, T>, size_t> b)
				  {return a.first.get(dim) > b.first.get(dim);}); // descending order!
		
		T minradius = all_eles.last().first.get(dim);
		//std::cout << "Min. radius: " << minradius << std::endl;
		
		max_level = std::ceil(std::log2(D_m / minradius)) - 1;
		// fun fact: more than about 64/dim levels are not possible with these indices (and thats already with 8 byte for a size_t). Perhaps set it as maximum? TODO.
		//std::cout << "Max. level: " << max_level << std::endl;
		
		auto end_iter = all_eles.begin();
		auto begin_iter(end_iter);
		unsigned int k = 0;
		
		// make entries in cells for the few particles that may live on the unpartitioned level 0 (in cell 0, yeah, mine are zero-based)
		while(end_iter->first.get(dim) >= D_m / (1l << (k+1)))
			++end_iter;
		if (begin_iter != end_iter)
			cells.insert(std::make_pair(0, std::make_pair(begin_iter, end_iter)));
		
		// now for the "normal" levels where we do the partition-based sorting:
		for (k = 1; k <= max_level; ++k) {
			// We are populating level k with the particles that belong here
			// Move iterators so that only these particles are in the range
			begin_iter = end_iter;
			while(end_iter != all_eles.end() && end_iter->first.get(dim) >= D_m / (1l << (k+1)))
				++end_iter;
			//std::cout << "Sorting on level " << k << ": " << std::distance(begin_iter, end_iter) << " of " << all_eles.size() << std:: endl;
			
			// If there are no particles on this level... no need to fill the tree/hashmap.
			if (begin_iter == end_iter)
				continue;
			
			// Spatial sorting
			// We can tolerate the recursion, it only goes k levels down.
			sort(begin_iter, end_iter, 1, k, 0, sbox, 0, 0);
		}
	}
	
	std::pair<std::string, std::string> printTree(unsigned int current_level, size_t parent_row_no, size_t current_sibling_no) {
		bool foundsomething = false;
		std::stringstream result;
		size_t index = getIndex(current_level, parent_row_no, current_sibling_no);
		auto resultiter = cells.find(index);
		for(unsigned int i=0; i<current_level; i++)
			result << "  ";
		result << "(" << current_sibling_no << ") " << index;
		if(resultiter != cells.end()) {
			foundsomething = true;
			result << ": " << std::distance(resultiter->second.first, resultiter->second.second) << ": ";
			for(auto childiter = resultiter->second.first; childiter != resultiter->second.second; ++childiter)
				result << "(" << childiter->first.toString() << ") ";
		}
		result << std::endl;
		
		std::string childresult;
		if(current_level <= max_level)
			for(int i=0; i<(1l << dim); i++) {
				auto child = printTree(current_level+1, (parent_row_no << dim) + current_sibling_no, i);
				if(child.second != "")
					foundsomething = true;
				childresult += child.first + child.second;
			}
		
		return std::make_pair(result.str(), foundsomething ? childresult : "");
	}
	
	/*! \brief called Op1 in paper
	 */
	inline size_t findCellIndex(Point<dim+1, T> p)
	{
		T radius = p.get(dim);
		if(radius > D_m)
			throw std::invalid_argument("Point with cutoff radius > D_m searched!");
			// Alternatively we could just limit it to D_m, but... that might yield wrong results.
		
		unsigned int target_level = std::ceil(std::log2(D_m / radius)) - 1;
		//std::cout << p.toString() << " is on level " << target_level << std::endl;
		
		if(target_level == 0)
			return 0;
		
		size_t cellindexoffsetinlastlevel = 0;
		SpaceBox<dim,T> cellbox(sbox);
		
		for (unsigned int k = 1; k <= target_level; ++k) {
			// Which part of the partition are we in? (0 - 2^dim)
			for (unsigned int i = 0; i < dim; ++i) {
				const T pivot = (cellbox.getHigh(i) + cellbox.getLow(i)) / 2;
				//std::cout << "Check in dim " << i << " coord " << p.get(i) << " against " << pivot << std::endl;
				if(p.get(i) >= pivot) {
					cellindexoffsetinlastlevel += (1l << static_cast<size_t>((dim-1 -i) + dim * (target_level - k)));
					// first summand gives shift in current level, second summand multiplies that into last level (sum in exponent is multiplication):
					// cellindexoffsetinlastlevel += 2 ^ (dim-1 - i) * 2 ^ dim * ... * 2 ^ dim -- (target_level - k) times 
					cellbox.setLow(i, pivot);
				} else
					cellbox.setHigh(i, pivot);
			}
			
			//std::cout << "In k=" << k << ": firstindex=" << firstindexoflastlevel << " offset=" << cellindexoffsetinlastlevel << std::endl;
		}
		
		return firstIndexOfLevel(target_level) + cellindexoffsetinlastlevel;
	}
	
	/*! \brief gets start and end iterator of contents of a cell (or any two identical iterators, if the cell is empty)
	 */
	inline std::pair<decltype(all_eles.begin()), decltype(all_eles.end())> getCellContents(size_t cellindex)
	{
		auto iter = cells.find(cellindex);
		if (iter == cells.end())
			return std::make_pair(all_eles.end(), all_eles.end()); // nothing in there
		else
			return std::make_pair(iter->second.first, iter->second.second);
	}
	
	/*! \brief called Op3 in paper
	 */
	inline std::vector<size_t> findChildrenIndices(size_t rootindex)
	{
		/*
		 * ((J+1) * 2^d + l) - 1 = J ^ 2^d + (l + 2^d - 1) -- J+1 and the global -1 because our indices start at 0 and stuff.
		 * 
		 * so l is not in [-(2^d)+2 ; 1], but in [1 ; 2^d]
		 */
		std::vector<size_t> children;
		children.reserve(1 << dim);
		for(size_t l = 1; l <= (1l << dim); ++l)
			children.push_back((rootindex << dim) + l);
		return children;
	}
	
	/*! \brief needed for Op5
	 */
	inline Point<dim+1,T> findCellCenter(size_t cellindex)
	{
		// inverse of (sum 2^(dim*i) over i from 0 to (level-1))  ->  level = floor(ld(index * (2^dim - 1) + 1) / dim)
		size_t level = std::floor(std::log2(((1l << dim) - 1l) * cellindex + 1l) / dim);
		
		size_t indexoffset = cellindex - firstIndexOfLevel(level); // offset in last level
		size_t indexmax = 1l << (dim*level); // also in last level
		
		//std::cout << "Index " << cellindex << " has offset " << indexoffset << "/" << indexmax << " in level " << level << std::endl;
		
		Point<dim+1,T> pos;
		// spatial coords
		for(unsigned int d = 0; d < dim; d++)
			pos.get(d) = sbox.getP1().get(d);
		// cutoff radius
		pos.get(dim) = (D_m / static_cast<T>(1l << (level+1))) * static_cast<T>(1.5f); // the *1.5 should be useless. TODO: make sure that it is
		
		// using this offset as a pivot on the whole row gives me the position w.r.t to the highest dimension on the highest level
		// then pivoting on this half of the row gives me the pos in the next highest dim and so on until
		// pivoting gives me the position in the highest dimension on the second highest level... and so on!
		
		size_t pivot = indexmax / 2l;
		size_t iter = 1;
		while(pivot > 0) {
			for(unsigned int d = 0; d < dim && pivot > 0; ++d) {
				//std::cout << "Pivot " << pivot << std::endl;
				if (indexoffset >= pivot) {
					indexoffset -= pivot;
					pos.get(d) += (sbox.getHigh(d) - sbox.getLow(d)) / (1 << iter);
				}
				pivot /= 2l;
			}
			++iter;
		}
		
		for(unsigned int d = 0; d < dim; ++d)
			pos.get(d) += (sbox.getHigh(d) - sbox.getLow(d)) / (1l << (level+1));
		
		return pos;
	}
	
	/*! \brief Get an element in the cell
	 *
	 * \param ele_id element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t& get(size_t ele_id)
	{
		return all_eles.get(ele_id).second;
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	inline void swap(AdaptiveCellList<dim,T,BALANCED,ele_container> & cl)
	{
		all_eles.swap(cl.all_eles);
	}

	/*! \brief Get the Nearest Neighborhood iterator
	 *
	 * \param cell cell id
	 *
	 */
	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,FULL,impl> getNNIterator(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,FULL,impl> cln(*this);

		return cln;
	}

	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,SYM,impl> cln(*this);

		return cln;
	}

	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,ele_container>,CRS,impl> cln(*this);

		return cln;
	}
};


#endif /* CELLLISTSTANDARD_HPP_ */
