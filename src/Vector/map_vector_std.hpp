/*
 * map_vector_std.hpp
 *
 *  Created on: Mar 8, 2015
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_HPP_
#define MAP_VECTOR_STD_HPP_

#include "se_vector.hpp"

/*! \brief Implementation of 1-D std::vector like structure
 *
 * this implementation is just a wrapper for the std::vector in the case
 * of data where the members cannot be parsed see openFPM_data wiki for more information
 *
 * ### Create add and access the elements
 * \snippet vector_test_util.hpp Create add and access stl
 *
 * \param T base type
 *
 */
template<typename T>
class vector<T,HeapMemory,grow_policy_double,STD_VECTOR>
{
	//! Actual size of the vector, warning: it is not the space allocated in grid
	//! grid size increase by a fixed amount every time we need a vector bigger than
	//! the actually allocated space
	size_t v_size;

	//! 1-D static grid
	std::vector<T> base;

	//! Error code
	size_t err_code;

public:

	//! it define that it is a vector
	typedef int yes_i_am_vector;

	//! iterator for the vector
	typedef vector_key_iterator iterator_key;
	//! Type of the value the vector is storing
	typedef T value_type;

	//! return the size of the vector
	inline size_t size() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return base.size();
	}


	/*! \ brief Resize the vector to contain n elements
	 *
	 * \param slot number of elements
	 *
	 */
	inline void resize(size_t slot)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		v_size = slot;

		base.resize(slot);
	}

	/*! \brief Remove all the element from the vector
	 *
	 */
	inline void clear()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.clear();
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * \param v element to add
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *vector_isel<T>::value
	 */
	inline void add(const T & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.push_back(v);
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * \param v element to add
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *vector_isel<T>::value
	 */
	inline void add(T && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.emplace_back(v);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.emplace_back(T());
	}

	/*! \brief Erase the elements from start to end
	 *
	 * \param start element
	 * \param end element
	 *
	 */
	void erase(typename std::vector<T>::iterator start, typename std::vector<T>::iterator end)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.erase(start,end);
	}

	/*! \brief Remove one entry from the vector
	 *
	 * \param key element to remove
	 *
	 */
	void remove(size_t key)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		vector_overflow(key);
#endif
		base.erase(base.begin() + key);
	}

	/*! \brief Return an std compatible iterator to the first element
	 *
	 * \return an iterator to the first element
	 *
	 */
	inline auto begin() -> decltype(base.begin())
	{
		return base.begin();
	}

	/*! \brief Return an std compatible iterator to the last element
	 *
	 * \return an iterator to the last element
	 *
	 */
	inline auto end() -> decltype(base.begin())
	{
		return base.end();
	}

	/*! \brief Get the last element
	 *
	 * \return the last element as reference
	 *
	 */
	inline T & last()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		if (base.size() == 0)
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " vector of size 0\n";
#endif
		return base[base.size()-1];
	}

	/*! \brief Duplicate the vector
	 *
	 * \return the duplicated vector
	 *
	 */
	std::vector<T> duplicate()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return base;
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param v vector to swap
	 *
	 */
	void swap(std::vector<T> && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.swap(v);
	}

	/*! \brief It eliminate double entries
	 *
	 * \note The base object must have an operator== defined
	 *
	 */
	void unique()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		auto it = std::unique(base.begin(),base.end());
		base.resize( std::distance(base.begin(),it) );
	}

	/*! \brief It sort the vector
	 *
	 * \note The base object must have an operator< defined
	 *
	 */
	void sort()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		std::sort(base.begin(), base.end());
	}

	/*! \brief Get an element of the vector
	 *
	 * \tparam p must be 0
	 *
	 * \param id element to get
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline T& get(size_t id)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		if (p != 0)
		{std::cerr << "Error the property does not exist" << "\n";}

		vector_overflow(id);
#endif

		return base[id];
	}

	/*! \brief Get an element of the vector
	 *
	 * \tparam p must be 0
	 *
	 * \param id element to get
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline const T& get(size_t id) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		if (p != 0)
		{std::cerr << "Error the property does not exist" << "\n";}

		vector_overflow(id);
#endif

		return base[id];
	}

	/*! \brief Get an element of the vector
	 *
	 * \param id element to get
	 *
	 * \return the element reference
	 *
	 */
	inline T & get(size_t id)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		if (id >= base.size())
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";
#endif
		return base[id];
	}

	/*! \brief Get an element of the vector
	 *
	 * \param id element to get
	 *
	 * \return the element value
	 *
	 */
	inline const T & get(size_t id) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		vector_overflow(id);
#endif
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
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		memset(&base[0],fl,base.size() * sizeof(T));
	}

	/*! \brief reserve a memory space in advance to avoid reallocation
	 *
	 * \param ns number of element the memory has to store
	 *
	 */

	inline void reserve(size_t ns)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.reserve(ns);
	}

	//! Constructor, vector of size 0
	vector() noexcept
	:v_size(0),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor, vector of size sz
	vector(size_t sz) noexcept
	:v_size(0),base(sz),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor from another vector
	vector(const vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & v) noexcept
	:v_size(0),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif

		base = v.base;
	}

	//! Constructor from another vector
	vector(vector<T,HeapMemory,grow_policy_double,STD_VECTOR> && v) noexcept
	:v_size(0),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif

		base.swap(v.base);
	}

	//! destructor
	~vector() noexcept
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.swap(v.base);
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & operator=(const vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base = v.base;

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & operator=(vector<T,HeapMemory,grow_policy_double,STD_VECTOR> && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.swap(v.base);

		return *this;
	}

	/*! \brief Check that two vectors are equal
	 *
	 * \param vector to compare
	 *
	 */
	bool operator!=(const vector<T, HeapMemory,grow_policy_double,STD_VECTOR> & v) const
	{
		return base != v.base;
	}

	/*! \brief Check that two vectors are not equal
	 *
	 * \param vector to compare
	 *
	 */
	bool operator==(const vector<T, HeapMemory,grow_policy_double,STD_VECTOR> & v) const
	{
		return base == v.base;
	}

	/*! \brief Get iterator
	 *
	 * \return an iterator
	 *
	 */
	vector_key_iterator getIterator() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
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

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	void * getPointer()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return &base[0];
	}

	/*! \brief This class has pointer inside
	 *
	 * \return false
	 *
	 */
	static bool noPointers()
	{
		return false;
	}

	/*! \brief Return the last error
	 *
	 */
	size_t getLastError()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return err_code;
	}

	/*! \brief check that the id does not overflow the buffer
	 *
	 * \param id to check
	 *
	 */
	inline void vector_overflow(size_t v1) const
	{
		if (v1 >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << v1 << "\n";\
			size_t * err_code_pointer = (size_t *)&this->err_code;\
			*err_code_pointer = 2001;\
			ACTION_ON_ERROR(VECTOR_ERROR);\
		}
	}

	/* \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 */
	long int who()
	{
#ifdef SE_CLASS2
		return check_whoami(this,8);
#endif
	}
};


#endif /* MAP_VECTOR_STD_HPP_ */
