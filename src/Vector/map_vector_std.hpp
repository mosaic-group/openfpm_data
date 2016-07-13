/*
 * map_vector_std.hpp
 *
 *  Created on: Mar 8, 2015
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_HPP_
#define MAP_VECTOR_STD_HPP_

#include "se_vector.hpp"
#include "map_vector_std_ptr.hpp"

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
class vector<T,HeapMemory,typename memory_traits_lin<T>::type,memory_traits_lin,grow_policy_double,STD_VECTOR>
{
	// Memory layout
	typedef typename memory_traits_lin<T>::type layout;
//      Doeas not work on gcc 4.8.4
//	template <typename lb> using layout_base = memory_traits_lin<lb>;

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

	typedef void base_to_copy;

	//This file implements a pack and unpack for std vector
#include "vector_std_pack_unpack.ipp"

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

		// here we have to check if the vector go into reallocation
		void * ptr_old = &base[0];
#endif

		base.resize(slot);

#ifdef SE_CLASS2
		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],slot*sizeof(T),VECTOR_STD_EVENT,1);
		}
#endif
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
		void * ptr_old = &base[0];
#endif

		base.push_back(v);

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
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
		void * ptr_old = &base[0];
#endif

		base.emplace_back(v);

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		base.emplace_back(T());

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
	}

	/*! \brief add elements to the vector
	 *
	 * \param eles elements to add
	 *
	 */
	template<typename Mem,typename l,template<typename> class lb,typename gp> inline void add(const openfpm::vector<T,Mem,l,lb,gp> & eles)
	{

#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		size_t start = base.size();
		base.resize(base.size() + eles.size());

		// copy the elements
		std::copy(eles.begin(),eles.end(),base.begin()+start);

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
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
	template<typename S> inline void add(const S & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		push_back_op<is_vector<T>::value,is_vector<S>::value,T,S>::push_back(base,v);

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
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
	template<typename S> inline void add(const S && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		base.push_back(v);

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
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

	/*! \brief Return an std compatible iterator to the first element
	 *
	 * \return an iterator to the first element
	 *
	 */
	inline auto begin() const -> const decltype(base.begin())
	{
		return base.begin();
	}

	/*! \brief Return an std compatible iterator to the last element
	 *
	 * \return an iterator to the last element
	 *
	 */
	inline auto end() const -> const decltype(base.begin())
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
	openfpm::vector<T> duplicate() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return *this;
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
		vector_overflow(id);
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
	:err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor, vector of size sz
	vector(size_t sz) noexcept
	:base(sz),err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
		check_new(&base[0],sizeof(T)*sz,VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor from another vector
	vector(const vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v) noexcept
	:err_code(0)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
		void * ptr_old = &base[0];
#endif

		base = v.base;

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif
	}

	/*! \brief Initializer from constructor
	 *
	 * \param v Initializer list
	 *
	 */
	vector(const std::initializer_list<T> & v)
	:base(v)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_STD_EVENT,1);
		check_new(&base[0],sizeof(T)*v.size(),VECTOR_STD_EVENT,1);
#endif
	}

	//! Constructor from another vector
	vector(vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v) noexcept
	:err_code(0)
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
		check_delete(&base[0]);
#endif
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.swap(v.base);
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v)
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
	vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(const vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		base = v.base;

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	template<typename Mem, typename gp> vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(const vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> & v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
		void * ptr_old = &base[0];
#endif

		base_copy<has_base_to_copy<vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR>>::value,
		          decltype(*this),
				  vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> >::copy(*this,v);
//		base = v.base;

#ifdef SE_CLASS2

		if (ptr_old != &base[0])
		{
			check_delete(ptr_old);
			check_new(&base[0],base.size()*sizeof(T),VECTOR_STD_EVENT,1);
		}

#endif

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		base.swap(v.base);

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	template<typename Mem, typename gp>  vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> && v)
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
	bool operator!=(const vector<T, HeapMemory, layout, memory_traits_lin,grow_policy_double,STD_VECTOR> & v) const
	{
		return base != v.base;
	}

	/*! \brief Check that two vectors are not equal
	 *
	 * \param vector to compare
	 *
	 */
	bool operator==(const vector<T, HeapMemory, layout, memory_traits_lin,grow_policy_double,STD_VECTOR> & v) const
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

	/*! \brief Get iterator untill a specified key
	 *
	 * \return an iterator
	 *
	 */
	vector_key_iterator getIteratorTo(size_t k) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return vector_key_iterator(k);
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
	template<int ... prp> inline size_t packMem(size_t n, size_t e) const
	{
		if (n == 0)
			return 0;
		else {
#ifdef DEBUG
			std::cout << "Inside map_vector_std.hpp packMem()" << std::endl;
#endif
			packMem_cond<has_packMem<T>::type::value, openfpm::vector<T, HeapMemory, layout, memory_traits_lin, grow_policy_double>, prp...> cm;
			return cm.packMemory(*this,n,0);
		}
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
	inline static size_t calculateMemDummy(size_t n, size_t e)
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

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	const void * getPointer() const
	{
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
#else
		return -1;
#endif
	}
};


#endif /* MAP_VECTOR_STD_HPP_ */
