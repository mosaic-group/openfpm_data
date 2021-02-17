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

#define OBJECT_ADD false
#define VECTOR_ADD true

//! struct to merge two vectors
template<bool objv, typename vect_dst>
struct add_prp_impl
{
	/*! \brief It add the element of a source vector to this vector
	 *
	 * The number of properties in the source vector must be smaller than the destination
	 * all the properties of S must be mapped so if S has 3 properties
	 * 3 numbers for args are required
	 *
	 * \tparam S Base object of the source vector
	 * \tparam M memory type of the source vector
	 * \tparam gp Grow policy of the source vector
	 * \tparam args one or more number that define which property to set-up
	 *
	 * \param v_src source vector
	 * \param v_dst destination vector
	 *
	 */
	template <typename S, typename M, typename gp, unsigned int impl, unsigned int ...args> inline static void add(const vector<S,M,typename memory_traits_lin<S>::type,memory_traits_lin,gp,impl> & v_src, vect_dst & v_dst)
	{
		//! Add the element of v
		for (size_t i = 0 ; i < v_src.size() ; i++)
		{
			// Add a new element
			v_dst.add();

			// equal object
			v_dst.get(v_dst.size()-1) = v_src.get(i);
		}
	}
};

//! struct to merge two vectors
template<typename vect_dst>
struct add_prp_impl<OBJECT_ADD,vect_dst>
{
	/*! \brief It add the element of a source vector to a destination vector
	 *
	 * The number of properties in the source vector must be smaller than the destination
	 * all the properties of S must be mapped so if S has 3 properties
	 * 3 numbers for args are required
	 *
	 * \tparam S Base object of the source vector
	 * \tparam M memory type of the source vector
	 * \tparam gp Grow policy of the source vector
	 * \tparam args one or more number that define which property to set-up
	 *
	 * \param v_src vector to merge
	 * \param v_dst vector to merge and result of the merge
	 *
	 */
	template <typename S, typename M, typename gp, unsigned int impl, unsigned int ...args> inline static void add(const vector<S,M,typename memory_traits_lin<S>::type,memory_traits_lin,gp,impl> & v_src, vect_dst & v_dst)
	{
			// Add a new element
			v_dst.add();

			// equal object
			v_dst.get(v_dst.size()-1) = v_src;
	}
};


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
template<typename T, typename grow_p>
class vector<T,HeapMemory,typename memory_traits_lin<T>::type,memory_traits_lin,grow_p,STD_VECTOR>
{
	// Memory layout
	typedef typename memory_traits_lin<T>::type layout;
//      Doeas not work on gcc 4.8.4
//	template <typename lb> using layout_base = memory_traits_lin<lb>;

	//! 1-D static grid
	std::vector<T> base;

	//! Error code
	size_t err_code = 0;

public:

	//! it define that it is a vector
	typedef int yes_i_am_vector;

	typedef memory_traits_lin<T> layout_base_;

	//! iterator for the vector
	typedef vector_key_iterator iterator_key;
	//! Type of the value the vector is storing
	typedef T value_type;

	typedef void base_to_copy;

	//! growing policy of this vector
	typedef grow_policy_double grow_policy;

	template<typename Tobj>
	struct layout_base__
	{
		typedef memory_traits_lin<Tobj> type;
	};

	//This file implements a pack and unpack for std vector
#include "vector_std_pack_unpack.ipp"

	//! return the size of the vector
	inline size_t size() const
	{
		return base.size();
	}


	/*! \ brief Resize the vector to contain n elements
	 *
	 * \param slot number of elements
	 *
	 */
	inline void resize(size_t slot)
	{
		base.resize(slot);
	}

	/*! \brief Remove all the element from the vector
	 *
	 */
	inline void clear()
	{
		base.clear();
	}

	/*! \brief Fit the memory to the size of the vector
	 *
	 */
	inline void shrink_to_fit()
	{
		base.shrink_to_fit();
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
		if (std::is_same<grow_p,openfpm::grow_policy_identity>::value == true)
		{
			// we reserve just one space more to avoid the capacity to increase by two
			base.reserve(base.size()+1);
		}

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
		if (std::is_same<grow_p,openfpm::grow_policy_identity>::value == true)
		{
			// we reserve just one space more to avoid the capacity to increase by two
			base.reserve(base.size()+1);
		}

		base.emplace_back(v);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */
	inline void add()
	{
		base.emplace_back(T());
	}

	/*! \brief add elements to the vector
	 *
	 * \param eles elements to add
	 *
	 */
	template<typename Mem,typename l,template<typename> class lb,typename gp> inline void add(const openfpm::vector<T,Mem,l,lb,gp> & eles)
	{
		if (std::is_same<grow_p,openfpm::grow_policy_identity>::value == true)
		{
			// we reserve just one space more to avoid the capacity to increase by two
			base.reserve(base.size() + eles.size());
		}

		size_t start = base.size();
		base.resize(base.size() + eles.size());

		// copy the elements
		std::copy(eles.begin(),eles.end(),base.begin()+start);
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the object
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
		push_back_op<is_vector<T>::value,is_vector<S>::value,T,S>::push_back(base,v);
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
		if (std::is_same<grow_p,openfpm::grow_policy_identity>::value == true)
		{
			// we reserve just one space more to avoid the capacity to increase by two
			base.reserve(base.size() + 1);
		}

		base.push_back(v);
	}

	/*! \brief It add the element of a source vector to this vector
	 *
	 * The number of properties in the source vector must be smaller than the destination
	 * all the properties of S must be mapped so if S has 3 properties
	 * 3 numbers for args are required
	 *
	 * \tparam S Base object of the source vector
	 * \tparam M memory type of the source vector
	 * \tparam gp Grow policy of the source vector
	 * \tparam args one or more number that define which property to set-up
	 *
	 * \param v source vector
	 *
	 */
	template <typename S,
	          typename M,
			  typename gp,
			  unsigned int impl,
			  template <typename> class layout_base,
			  unsigned int ...args>
	void add_prp(const vector<S,M,typename layout_base<S>::type,layout_base,gp,impl> & v)
	{
		add_prp_impl<std::is_same<S,T>::value,typename std::remove_pointer<decltype(*this)>::type>::template add<S,M,gp,impl,args...>(v,*this);
	}

	/*! \brief It add the element of a source vector to this vector
	 *
	 * The number of properties in the source vector must be smaller than the destination
	 * all the properties of S must be mapped so if S has 3 properties
	 * 3 numbers for args are required
	 *
	 * \tparam S Base object of the source vector
	 * \tparam M memory type of the source vector
	 * \tparam gp Grow policy of the source vector
	 * \tparam args one or more number that define which property to set-up
	 *
	 * \param v source vector
	 *
	 */
	template <typename S,
	          typename M,
			  typename gp,
			  unsigned int impl,
			  template <typename> class layout_base,
			  unsigned int ...args>
	void add_prp_device(const vector<S,M,typename layout_base<S>::type,layout_base,gp,impl> & v)
	{
		add_prp_impl<std::is_same<S,T>::value,typename std::remove_pointer<decltype(*this)>::type>::template add<S,M,gp,impl,args...>(v,*this);
	}

	/*! \brief It add the element it is equivalent to add
	 *
	 * exist to respect a general interface template parameters are unused the explanation
	 * refer to the interface specification, but is unused in this case
	 *
	 * \tparam S Base object of the source vector
	 * \tparam M memory type of the source vector
	 * \tparam gp Grow policy of the source vector
	 * \tparam args one or more number that define which property to set-up
	 *
	 * \param v source vector
	 *
	 */
	template <typename S,
	          typename M,
			  typename gp,
			  unsigned int impl,
			  template <typename> class layout_base,
			  unsigned int ...args>
	void add_prp(const T & v)
	{
		add(v);
	}

	/*! \brief Erase the elements from start to end
	 *
	 * \param start element
	 * \param end element
	 *
	 */
	void erase(typename std::vector<T>::iterator start, typename std::vector<T>::iterator end)
	{
		base.erase(start,end);
	}

	/*! \brief Remove one entry from the vector
	 *
	 * \param key element to remove
	 *
	 */
	void remove(size_t key)
	{
#ifdef SE_CLASS1
		vector_overflow(key);
#endif
		base.erase(base.begin() + key);
	}

	/*! \brief Remove several entries from the vector
	 *
	 * \warning the keys in the vector MUST be sorted
	 *
	 * \param keys objects id to remove
	 * \param start key starting point
	 *
	 */
	void remove(openfpm::vector<size_t> & keys, size_t start = 0)
	{
		// Nothing to remove return
		if (keys.size() <= start )
			return;

		size_t a_key = start;
		size_t d_k = keys.get(a_key);
		size_t s_k = keys.get(a_key) + 1;

		// keys
		while (s_k < size())
		{
			// s_k should always point to a key that is not going to be deleted
			while (a_key+1 < keys.size() && s_k == keys.get(a_key+1))
			{
				a_key++;
				s_k = keys.get(a_key) + 1;
			}

			// In case of overflow
			if (s_k >= size())
				break;

			base[d_k] = base[s_k];
			d_k++;
			s_k++;
		}

		// re-calculate the vector size

		base.resize(base.size() - (keys.size() - start));
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
		return *this;
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param v vector to swap
	 *
	 */
	void swap(std::vector<T> && v)
	{
		base.swap(v);
	}

	/*! \brief It eliminate double entries
	 *
	 * \note The base object must have an operator== defined
	 *
	 */
	void unique()
	{
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
		memset(&base[0],fl,base.size() * sizeof(T));
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
	vector() noexcept
	:err_code(0)
	{
	}

	//! Constructor, vector of size sz
	vector(size_t sz) noexcept
	:base(sz),err_code(0)
	{
	}

	//! Constructor from another vector
	vector(const vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v) noexcept
	:err_code(0)
	{
		base = v.base;
	}

	/*! \brief Initializer from constructor
	 *
	 * \param v Initializer list
	 *
	 */
	vector(const std::initializer_list<T> & v)
	:base(v)
	{
	}

	//! Constructor from another vector
	vector(vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v) noexcept
	:err_code(0)
	{
		base.swap(v.base);
	}

	//! destructor
	~vector() noexcept
	{
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v)
	{
		base.swap(v.base);
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v)
	{
		base.swap(v.base);
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(const vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & v)
	{
		base = v.base;

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \return itself
	 *
	 */
	template<typename Mem, typename gp> vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(const vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> & v)
	{
		base_copy<has_base_to_copy<vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR>>::value,
		          decltype(*this),
				  vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> >::copy(*this,v);

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> && v)
	{
		base.swap(v.base);

		return *this;
	}

	/*! \brief Operator= copy the vector into another
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	template<typename Mem, typename gp>  vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR> & operator=(vector<T,Mem,layout,memory_traits_lin,gp,STD_VECTOR> && v)
	{
		base.swap(v.base);

		return *this;
	}

	/*! \brief Check that two vectors are equal
	 *
	 * \param v vector to compare
	 *
	 * \return true if they differs
	 *
	 */
	bool operator!=(const vector<T, HeapMemory, layout, memory_traits_lin,grow_policy_double,STD_VECTOR> & v) const
	{
		return base != v.base;
	}

	/*! \brief Check that two vectors are not equal
	 *
	 * \param v vector to compare
	 *
	 * \return true if the vector match
	 *
	 */
	bool operator==(const vector<T, HeapMemory, layout, memory_traits_lin,grow_policy_double,STD_VECTOR> & v) const
	{
		return base == v.base;
	}

	/*! \brief Get an iterator over all the elements of the vector
	 *
	 * \return an iterator
	 *
	 */
	vector_key_iterator getIterator() const
	{
		return vector_key_iterator(base.size());
	}

	/*! \brief Get iterator until a specified element
	 *
	 * \param k key
	 *
	 * \return an iterator
	 *
	 */
	vector_key_iterator getIteratorTo(size_t k) const
	{
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
		else
		{
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
		return &base[0];
	}

	/*! \brief Do nothing
	 *
	 */
	template<unsigned int ... prp> void hostToDevice()
	{}

	/*! \brief Do nothing
	 *
	 */
	template<unsigned int ... prp> void deviceToHost()
	{}

	/*! \brief Do nothing
	 *
	 *
	 */
	template<unsigned int ... prp> void deviceToHost(size_t start, size_t stop)
	{}

	/*! \brief Do nothing
	 *
	 *
	 */
	template<unsigned int ... prp> void hostToDevice(size_t start, size_t stop)
	{}

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
	 * \return erro code
	 *
	 */
	size_t getLastError()
	{
		return err_code;
	}

	/*! \brief check that the id does not overflow the buffer
	 *
	 * \param v1 id to check
	 *
	 */
	inline void vector_overflow(size_t v1) const
	{
		if (v1 >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << v1 << "\n";\
			size_t * err_code_pointer = (size_t *)&this->err_code;\
			*err_code_pointer = 2001;\

			ACTION_ON_ERROR(VECTOR_ERROR_OBJECT);\
		}
	}
};

/*! \brief Implementation of 1-D std::vector like structure
 *
 * this implementation is just a wrapper for the std::vector. It work a little different from vector.
 * In general for a normal vector of objects A vector<A> if you resize to zero, the destructor  of
 * the object A is called.This vector differ in this behaviour. the destructor is not called. This give the possibility
 * to have a set of fully retained objects. This class is just a simple wrapper for the normal openfpm::vector where
 * size and resize are redefined to change the behaviour. A destructive resize is callable with resize_base(), and the internal
 * size of the base vactor can be queried with size_base()
 *
 * \param T base type
 *
 */
template<typename T>
class vector_fr
:private vector<T,HeapMemory,typename memory_traits_lin<T>::type,memory_traits_lin,grow_policy_double,STD_VECTOR>
{
	typedef vector<T,HeapMemory,typename memory_traits_lin<T>::type,memory_traits_lin,grow_policy_double,STD_VECTOR> base_type;

	//! size of the vector
	size_t v_size = 0;

public:

	/*! \brief return the size of the vector
	 *
	 * \return the size
	 *
	 */
	size_t size()
	{
		return v_size;
	}

	/*! \brief return the base size of the vector
	 *
	 * \return the size
	 *
	 */
	size_t size_base()
	{
		return base_type::size();
	}

	/*! \brief resize the vector retaining the objects
	 *
	 * \param new size of the vector
	 *
	 */
	void resize(size_t sz)
	{
		if (sz <= base_type::size())
		{
			v_size = sz;
			return;
		}

		base_type::resize(sz);

		v_size = sz;
	}

	/*! \brief Eliminate all elements
	 *
	 *
	 */
	void clear()
	{
		resize(0);
	}

	/*! \brief Add another element
	 *
	 *
	 */
	void add()
	{
		resize(v_size+1);
	}

	/*! \brief resize the base vector (this kill the objects)
	 *
	 * \param new size of the vector
	 *
	 */
	void resize_base(size_t sz)
	{
		base_type::resize(sz);
		v_size = sz;
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
		return base_type::get(id);
	}

	/*! \brief Get an element of the vector
	 *
	 * \param id element to get
	 *
	 * \return the element reference
	 *
	 */
	inline T & last()
	{
		return base_type::get(size()-1);
	}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector_fr<T> & v)
	{
		size_t v_size_tmp = v.v_size;
		v.v_size = v_size;
		v_size = v_size_tmp;

		base_type::swap(v);
	}
};

#endif /* MAP_VECTOR_STD_HPP_ */
