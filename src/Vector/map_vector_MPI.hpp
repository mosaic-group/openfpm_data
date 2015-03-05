/*
 * map_vector_MPI_Request.hpp
 *
 *  Created on: Mar 4, 2015
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_MPI_REQUEST_HPP_
#define MAP_VECTOR_MPI_REQUEST_HPP_

#ifdef HAVE_MPI

/*! \brief specialization for MPI_Request
 *
 */

template<>
struct device_cpu<MPI_Request>
{};

/*! \brief specialization for MPI_Status
 *
 */

template<>
struct device_cpu<MPI_Status>
{};

/*! \brief Implementation of 1-D std::vector like structure
 *
 * this implementation is just a wrapper for the std::vector in the case
 * of the primitive MPI_Request
 *
 * \param T base type
 *
 */

template<typename Memory>
class vector<MPI_Request,device_cpu<MPI_Request>,Memory>
{
	//! 1-D static grid
	std::vector<MPI_Request> base;

public:

	// iterator for the vector
	typedef vector_key_iterator iterator_key;

	//! return the size of the vector
	inline size_t size()
	{
		return base.size();
	}


	/*! \brief Resize the vector
	 *
	 * Resize the vector
	 *
	 * \param how many slot to reserve
	 *
	 */

	inline void resize(size_t slot)
	{
		base.resize(slot);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
		base.resize(base.size() + 1);
	}

	/*! \brief Get the last element
	 *
	 * \return the last element
	 *
	 */
	inline MPI_Request & last()
	{
		return base[base.size()-1];
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *
	 */
	inline void add(const MPI_Request & v)
	{
		base.push_back(v);
	}

	/*! \brief Duplicate the vector
	 *
	 * \return the duplicated vector
	 *
	 */
	std::vector<MPI_Request> duplicate()
	{
		return base;
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param vector to swap
	 *
	 */
	void swap(std::vector<MPI_Request> && v)
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
	template <unsigned int p>inline MPI_Request & get(size_t id)
	{
	#ifdef DEBUG
		if (p != 0)
		{std::cerr << "Error the property does not exist" << "\n";}

		if (id >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << "  " << __LINE__ << " overflow id: " << id << "\n";
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
	inline MPI_Request & get(size_t id)
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

	void swap(openfpm::vector<MPI_Request,device_cpu<MPI_Request>,Memory> & v)
	{
		base.swap(v.base);
	}

	/*! \brief Get the element i
	 *
	 * \param i element to get
	 * \return the reference to the element i
	 *
	 */
	MPI_Request & operator[](size_t i)
	{
		return base[i];
	}
};

template<typename Memory>
class vector<MPI_Status,device_cpu<MPI_Status>,Memory>
{
	//! 1-D static grid
	std::vector<MPI_Status> base;

public:

	// iterator for the vector
	typedef vector_key_iterator iterator_key;

	//! return the size of the vector
	inline size_t size()
	{
		return base.size();
	}


	/*! \brief Resize the vector
	 *
	 * Resize the vector
	 *
	 * \param how many slot to reserve
	 *
	 */

	inline void resize(size_t slot)
	{
		base.resize(slot);
	}

	/*! \brief Add an empty object (it call the default constructor () ) at the end of the vector
	 *
	 */

	inline void add()
	{
		base.resize(base.size() + 1);
	}

	/*! \brief Get the last element
	 *
	 * \return the last element
	 *
	 */
	inline MPI_Status & last()
	{
		return base[base.size()-1];
	}

	/*! \brief It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * It insert a new object on the vector, eventually it reallocate the grid
	 *
	 * \warning It is not thread safe should not be used in multi-thread environment
	 *          reallocation, work only on cpu
	 *
	 *
	 */
	inline void add(const MPI_Status & v)
	{
		base.push_back(v);
	}

	/*! \brief Duplicate the vector
	 *
	 * \return the duplicated vector
	 *
	 */
	std::vector<MPI_Status> duplicate()
	{
		return base;
	}

	/*! \brief swap the memory between the two vector
	 *
	 * \param vector to swap
	 *
	 */
	void swap(std::vector<MPI_Status> && v)
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
	template <unsigned int p>inline MPI_Status & get(size_t id)
	{
	#ifdef DEBUG
		if (p != 0)
		{std::cerr << "Error the property does not exist" << "\n";}

		if (id >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << "  " << __LINE__ << " overflow id: " << id << "\n";
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
	inline MPI_Status & get(size_t id)
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

	void swap(openfpm::vector<MPI_Status,device_cpu<MPI_Status>,Memory> & v)
	{
		base.swap(v.base);
	}

	/*! \brief Get the element i
	 *
	 * \param i element to ge
	 * \return the reference to the element i
	 *
	 */
	MPI_Status & operator[](size_t i)
	{
		return base[i];
	}
};

#endif

#endif /* MAP_VECTOR_MPI_REQUEST_HPP_ */
