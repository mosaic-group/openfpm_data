#ifndef ADAPTIVECELLLIST_HPP_
#define ADAPTIVECELLLIST_HPP_

#include "../CellList/CellList.hpp"
#include "../CellList/CellDecomposer.hpp"
#include "AdaptiveCellListNNIterator.hpp"
#include "Space/SpaceBox.hpp"


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
template<unsigned int dim, typename T,  unsigned int impl=BALANCED, typename ele_container=openfpm::vector<size_t>>
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
class AdaptiveCellList<dim,T,BALANCED,ele_container> // : public CellDecomposer_sm<dim,T>
{
	// vector of all elements
	ele_container all_eles;

	//A point
	Point<dim,T> orig;

        
        
public:

	// Object type that the structure store
	typedef T value_type;

        inline size_t size() {return all_eles.size();}
        
	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	void Initialize(SpaceBox<dim,T> & sbox, Point<dim,T> & orig)
	{
		this->orig = orig;
	}

	/*! \brief Cell list
	 *
	 * \param ... (Non default constructor if needed)
	 *
	 */
	AdaptiveCellList(SpaceBox<dim,T> & box, Point<dim,T> & orig)
	{
		Initialize(box,orig);
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate +1 (the last is the radius of the particle)
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim+1], typename ele_container::value_type ele)
	{
                all_eles.add(ele);
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim+1,T> & pos, typename ele_container::value_type ele)
	{
                all_eles.add(ele);
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
