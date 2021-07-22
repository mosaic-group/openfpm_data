/*
 * map_vector_std_cuda.hpp
 *
 *  Created on: Mar 7, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_CUDA_HPP_
#define MAP_VECTOR_STD_CUDA_HPP_

#include "map_vector_std_cuda_ker.cuh"

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
class vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR>
{
	// Memory layout
	typedef typename memory_traits_inte<aggregate<T>>::type layout;

	//! 1-D static grid
	vector<aggregate<T>,CudaMemory,memory_traits_inte,grow_policy_double> base;

	//! Error code
	size_t err_code = 0;

public:

	//! it define that it is a vector
	typedef int yes_i_am_vector;

	typedef memory_traits_inte<T> layout_base_;

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
		typedef memory_traits_inte<Tobj> type;
	};

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
		base.resize_no_device(slot);
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
	 * \param v element to add
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *vector_isel<T>::value
	 */
	inline void add(const T & v)
	{
		base.add_no_device();
		base.template get<0>(size()-1) = v;
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
		base.add_no_device();
		base.template get<0>(size()-1).swap(v);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */
	inline void add()
	{
		base.add_no_device();
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
		return base.template get<0>(size()-1);
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
		return base.template get<0>(size()-1);
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param v vector to swap
	 *
	 */
	void swap(vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> && v)
	{
		base.swap(v.base);
	}

	inline T& at(int i)
	{
		return get(i);
	}

	/*! \brief Get an element of the vector
	 *
	 * \param id element to get
	 *
	 * \return the reference to the element
	 *
	 */
	inline T& operator[](int id)
	{
#ifdef SE_CLASS1
		vector_overflow(id);
#endif

		return base.template get<0>(id);
	}

	/*! \brief Get an element of the vector
	 *
	 * \param id element to get
	 *
	 * \return the reference to the element
	 *
	 */
	inline T& get(int id)
	{
#ifdef SE_CLASS1
		vector_overflow(id);
#endif

		return base.template get<0>(id);
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
	inline const T& get(int id) const
	{
#ifdef SE_CLASS1
		vector_overflow(id);
#endif

		return base.template get<0>(id);
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
		return base.template get<0>(id);
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
		return base.template get<0>(id);
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
	vector(const vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> & v) noexcept
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
	// :base(v)
	{
		// causes error if simply passed to base
		for (T t: v)
			add(t);
	}

	//! Constructor from another vector
	vector(vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> && v) noexcept
	:err_code(0)
	{
		base.swap(v.base);
	}

	//! Constructor from iterators
	template< class InputIt >
	vector(InputIt first, InputIt last)
	{
		while(first != last) {
			add(*first);
			++first;
		}
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
	void swap(vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> & v)
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
	vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> &
	operator=(const vector<T,CudaMemory,memory_traits_inte,grow_policy_double,STD_VECTOR> & v)
	{
		base = v.base;

		return *this;
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

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	void * getPointer()
	{
		return &base.template get<0>(0);
	}

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	const void * getPointer() const
	{
		return &base.template get<0>(0);
	}

	/*! \brief Do nothing
	 *
	 */
	template<unsigned int ... prp> void hostToDevice()
	{
		base.template hostToDevice<0>();
	}

	/*! \brief Do nothing
	 *
	 */
	template<unsigned int ... prp> void deviceToHost()
	{
		base.template deviceToHost<0>();
	}

	/*! \brief Do nothing
	 *
	 *
	 */
	template<unsigned int ... prp> void deviceToHost(size_t start, size_t stop)
	{
		base.template deviceToHost<0>(start,stop);
	}

	/*! \brief Do nothing
	 *
	 *
	 */
	template<unsigned int ... prp> void hostToDevice(size_t start, size_t stop)
	{
		base.template hostToDevice<0>(start,stop);
	}

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 * \return an usable vector in the kernel
	 *
	 */
	vector_custd_ker<typename apply_transform<memory_traits_inte,aggregate<T>>::type,memory_traits_inte> toKernel()
	{
		vector_custd_ker<typename apply_transform<memory_traits_inte,aggregate<T>>::type,memory_traits_inte> v(base.size(),this->base.toKernel());

		return v;
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


#endif /* MAP_VECTOR_STD_CUDA_HPP_ */
