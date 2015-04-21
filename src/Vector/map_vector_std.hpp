/*
 * map_vector_std.hpp
 *
 *  Created on: Mar 8, 2015
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_HPP_
#define MAP_VECTOR_STD_HPP_


/*! \brief Implementation of 1-D std::vector like structure
 *
 * this implementation is just a wrapper for the std::vector in the case
 * of the primitive size_t
 *
 * \param T base type
 *
 */

template<typename T>
class vector<T,device_cpu<T>,HeapMemory,grow_policy_double,STD_VECTOR>
{
	//! Actual size of the vector, warning: it is not the space allocated in grid
	//! grid size increase by a fixed amount every time we need a vector bigger than
	//! the actually allocated space
	size_t v_size;

	//! 1-D static grid
	std::vector<T> base;

public:

	// iterator for the vector
	typedef vector_key_iterator iterator_key;
	//! Type of the value the vector is storing
	typedef T value_type;

	//! return the size of the vector
	inline size_t size()
	{
		return base.size();
	}


	/*! \ brief Resize the vector
	 *
	 * Resize the vector
	 *
	 * \param how many slot to reserve
	 *
	 */

	inline void resize(size_t slot)
	{
		v_size = slot;

		base.resize(slot);
	}

	/*! \brief Remove all the element from the vector
	 *
	 */
	inline void clear()
	{
		base.clear();
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *vector_isel<T>::value
	 */
	inline void add(const T & v)
	{
		base.push_back(v);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
		base.resize(base.size() + 1);
	}

	/*! \brief Remove one entry from the vector
	 *
	 * \param keys element to remove
	 *
	 */
	void remove(size_t key)
	{
		base.erase(base.begin() + key);
	}

	/*! \brief Get the last element
	 *
	 * \return the last element
	 *
	 */
	inline T & last()
	{
		return base[base.size()-1];
	}

	/*! \brief Duplicate the vector
	 *
	 * \return the duplicated vector
	 *
	 */

	std::vector<T> duplicate()
	{
		return base;
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param vector to swap
	 *
	 */

	void swap(std::vector<T> && v)
	{
		base.swap(v);
	}

	/*! \brief Get an element of the vector
	 *
	 * Get an element of the vector
	 *
	 * \tparam must be 0
	 *
	 * \param id Element to get
	 * \param p Property to get
	 *
	 */
	template <unsigned int p>inline size_t & get(size_t id)
	{
#ifdef DEBUG
		if (p != 0)
		{std::cerr << "Error the property does not exist" << "\n";}

		if (id >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";
		}
#endif

		return base[id];
	}

	/*! \brief Get an element of the vector
	 *
	 * Get an element of the vector
	 *
	 * \param id Element to get
	 * \param p Property to get
	 *
	 */
	inline T & get(size_t id)
	{
		return base[id];
	}

	/*! \brief it fill all the memory of fl patterns
	 *
	 * WARNING does not assign a value to each element but it fill the memory
	 * Useful to fast set the memory to zero
	 *
	 * \param fl byte to fill
	 *
	 */

	inline void fill(unsigned char fl)
	{
		memset(&base[0],0,base.size());
	}

	/*! \brief reserve a memory space in advance to avoid reallocation
	 *
	 * \param ns number of element the memory has to store
	 *
	 */

	inline void reserve(size_t ns)
	{
		base.reserve(ns);
	}

	//! Constructor, vector of size 0
	vector() {}

	//! Constructor, vector of size sz
	vector(size_t sz):base(sz) {}

	/*! swap the content of the vector
	 *
	 * \param vector to be swapped with
	 *
	 */

	void swap(openfpm::vector<T,device_cpu<T>,HeapMemory,grow_policy_double,STD_VECTOR> & v)
	{
		base.swap(v.base);
	}

	/*! \brief Get iterator
	 *
	 * \return an iterator
	 *
	 */

	vector_key_iterator getIterator()
	{
		return vector_key_iterator(base.size());
	}

	/*! \brief Calculate the memory size required to allocate n elements
	 *
	 * Calculate the total size required to store n-elements in a vector
	 *
	 * \param n number of elements
	 * \param e unused
	 *
	 * \return the size of the allocation number e
	 *
	 */
	inline static size_t calculateMem(size_t n, size_t e)
	{
		return n*sizeof(T);
	}

	/*! \brief How many allocation are required to create n-elements
	 *
	 * \param n number of elements
	 *
	 * \return the number of allocations
	 *
	 */
	inline static size_t calculateNMem(size_t n)
	{
		return 1;
	}
};


#endif /* MAP_VECTOR_STD_HPP_ */
