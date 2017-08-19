/*
 * map_vector_std_ptr.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_PTR_HPP_
#define OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_PTR_HPP_


template<typename T,typename gp>
class vector<T,PtrMemory,typename memory_traits_lin<T>::type,memory_traits_lin,gp,STD_VECTOR>
{
	//! Memory layout
	typedef typename memory_traits_lin<T>::type layout;

	//! function for memory layout
	template <typename lb> using layout_base = memory_traits_lin<lb>;

	//! Actual size of the vector, warning: it is not the space allocated in grid
	//! grid size increase by a fixed amount every time we need a vector bigger than
	//! the actually allocated space
	size_t v_size;

	//! 1-D static grid
	PtrMemory * mem;

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
		return v_size;
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
		// resize is valid only if v_size is 0 and it match the size of PtrMemory
		if (slot > mem->size()/sizeof(T))
			std::cerr << __FILE__ << ":" << __LINE__ << " error: this vector cannot be bigger than " << mem->size()/sizeof(T) << " elements\n";
		v_size = slot;
	}

	/*! \brief Remove all the element from the vector
	 *
	 */
	inline void clear()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		v_size = 0;
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
		std::cerr << __FILE__ << ":" << __LINE__ << " error: you cannot add a new element to this vector \n";
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
		std::cerr << __FILE__ << ":" << __LINE__ << " error: you cannot add new element to this vector \n";
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		v_size++;

		if (v_size > mem->size()/sizeof(T))
			std::cerr << __FILE__ << ":" << __LINE__ << " error: you cannot resize a PtrMemory vector over the size stored by PtrMemory";
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
		std::cerr << __FILE__ << ":" << __LINE__ << " error: you cannot erase element from this vector \n";
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
		std::cerr << __FILE__ << ":" << __LINE__ << " error: you cannot remove elements from this vector \n";
	}

	/*! \brief Return an std compatible iterator to the first element
	 *
	 * \return an iterator to the first element
	 *
	 */
	inline T * begin()
	{
		return mem->getPointer();
	}

	/*! \brief Return an std compatible iterator to the past the last element
	 *
	 * \return an iterator to the last element
	 *
	 */
	inline T * end()
	{
		return &((T *)mem->getPointer())[v_size];
	}

	/*! \brief Return an std compatible iterator to the first element
	 *
	 * \return an iterator to the first element
	 *
	 */
	inline const T * begin() const
	{
		return (T *)mem->getPointer();
	}

	/*! \brief Return an std compatible iterator to the past the last element
	 *
	 * \return an iterator to the last element
	 *
	 */
	inline const T * end() const
	{
		return &((T *)mem->getPointer())[v_size];
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
		if (v_size == 0)
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " vector of size 0\n";
#endif
		return ((T *)mem->getPointer())[v_size-1];
	}

	/*! \brief Get the last element
	 *
	 * \return the last element as reference
	 *
	 */
	inline const T & last() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		if (v_size == 0)
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " vector of size 0" << std::endl;
#endif
		return ((T *)mem->getPointer())[v_size-1];
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

		return ((T *)mem->getPointer())[id];
	}

	inline void setMemory(PtrMemory & mem)
	{
		this->mem = &mem;
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

		return ((T *)mem->getPointer())[id];
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
		if (id >= v_size)
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";
#endif
		return ((T *)mem->getPointer())[id];
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
		return ((T *)mem->getPointer())[id];
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
		memset(mem->getPointer(),fl,v_size * sizeof(T));
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
	}

	//! Constructor, vector of size 0
	vector() noexcept
	:v_size(0),mem(NULL),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor, vector of size sz
	vector(size_t sz) noexcept
	:v_size(sz),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor from another vector
	vector(const vector<T,PtrMemory,layout,layout_base,gp,STD_VECTOR> & v) noexcept
	:v_size(0),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif

		std::cerr << __FILE__ << ":" << __LINE__ << " error: copy constructor is not supported by this vector \n";
	}


	//! Constructor from another vector
	vector(vector<T,PtrMemory,layout,layout_base,gp,STD_VECTOR> && v) noexcept
	:v_size(0),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif

		std::cerr << __FILE__ << ":" << __LINE__ << " error: copy constructor is not supported by this vector \n";
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
	void swap(openfpm::vector<T,PtrMemory,layout,layout_base,gp,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		std::cerr << __FILE__ << ":" << __LINE__ << " error: swap is not supported by this vector \n";
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,layout,layout_base,grow_policy_double,STD_VECTOR> & operator=(const vector<T,HeapMemory,layout,layout_base,grow_policy_double,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		std::cerr << __FILE__ << ":" << __LINE__ << " error: operator= is not supported by this vector \n";

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,layout,layout_base,grow_policy_double,STD_VECTOR> & operator=(vector<T,HeapMemory,layout,layout_base,grow_policy_double,STD_VECTOR> && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		std::cerr << __FILE__ << ":" << __LINE__ << " error: operator= is not supported by this vector \n";

		return *this;
	}

	/*! \brief Check that two vectors are equal
	 *
	 * \param vector to compare
	 *
	 * \return true if the two vector match
	 *
	 */
	bool operator!=(const vector<T, HeapMemory, layout, layout_base,grow_policy_double,STD_VECTOR> & v) const
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error: operator!= is not supported by this vector \n";

		return false;
	}

	/*! \brief Check that two vectors are not equal
	 *
	 * \param vector to compare
	 *
	 * \return true if the vector match
	 *
	 */
	bool operator==(const vector<T, HeapMemory, layout, layout_base,grow_policy_double,STD_VECTOR> & v) const
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error: operator== is not supported by this vector \n";

		return false;
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
		return vector_key_iterator(v_size);
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
		return mem->getPointer();
	}

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	const void * getPointer() const
	{
		return mem->getPointer();
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
	 * \return the last error
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
	 * \param v1 to check
	 *
	 */
	inline void vector_overflow(size_t v1) const
	{
		if (v1 >= v_size)
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << v1 << "\n";\
			size_t * err_code_pointer = (size_t *)&this->err_code;\
			*err_code_pointer = 2001;\
			ACTION_ON_ERROR(VECTOR_ERROR_OBJECT);\
		}
	}

	/* \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 * \return the allocation id
	 *
	 */
	long int who()
	{
#ifdef SE_CLASS2
		return check_whoami(this,8);
#else
		return -1;
#endif
	}
};


#endif /* OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_PTR_HPP_ */
