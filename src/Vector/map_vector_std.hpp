/*
 * map_vector_std.hpp
 *
 *  Created on: Mar 8, 2015
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_HPP_
#define MAP_VECTOR_STD_HPP_

#include "Pack_stat.hpp"
#include "memory/ExtPreAlloc.hpp"

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

public:

	//! it define that it is a vector
	typedef int yes_i_am_vector;

	typedef typename grid_cpu<1,typename T::type>::memory_conf memory_conf;

	//! iterator for the vector
	typedef vector_key_iterator iterator_key;
	//! Type of the value the vector is storing
	typedef T value_type;

	// Structures to check if inside of the object there is/are simple or complex structure/structures

	template<bool cond, typename T1, typename Memory1, int ... prp>
	struct nested_pack_cond
	{
		void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
		{
			//std::cout << "No pack() function inside!" << std::endl;
		}

	};

	template<typename T1, typename Memory1, int ... prp>
	struct nested_pack_cond<true, T1, Memory1, prp...>
	{
		void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
		{
		   for (int i = 0; i < obj.size(); i++) {
			   obj.get(i).template pack<prp...>(mem, sts);
		   }
		}
	};


	/*! \brief pack a vector selecting the properties to pack
	 *
	 * \param mem preallocated memory where to pack the vector
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */


	template<int ... prp> void pack(ExtPreAlloc<HeapMemory> & mem, openfpm::vector<T> & obj, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		// if no properties should be packed return
		if (sizeof...(prp) == 0)
			return;

		nested_pack_cond<has_Pack<T>::type::value, T, HeapMemory, prp...> dc;
		dc.packing(mem, obj, sts);


		// Sending property object
		typedef openfpm::vector<T> vctr;
		typedef object<typename object_creator<typename T::value_type::type,prp...>::type> prp_object;

		typedef openfpm::vector<prp_object,ExtPreAlloc<HeapMemory>,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());
		auto obj_it = obj.getIterator();

		while (obj_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_src;
			// destination object type
			typedef encapc<1,prp_object,typename dtype::memory_conf > encap_dst;

			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,ENCAP,prp...>(obj.get(obj_it.get()),dest.get(obj_it.get()));

			++obj_it;
		}

		// Update statistic
		sts.incReq();

	}

	/*! \brief Insert an allocation request into the vector
	 *
	 * \param obj vector object to pack
	 * \param requests vector
	 *
	 */
	template<int ... prp> void packRequest(openfpm::vector<T> & obj, std::vector<size_t> & v)
	{
		//if (sizeof...(prp) == 0)
		//{
			openfpm::vector<T> vect;

			// Calculate the required memory for packing
			size_t alloc_ele = vect.calculateMem(obj.size(),0);

			v.push_back(alloc_ele);
		/*}
		else
		{
			typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;
			openfpm::vector<prp_object> vect;

			// Calculate the required memory for packing
			size_t alloc_ele = vect.calculateMem(obj.size(),0);

			v.push_back(alloc_ele);
		}*/
	}

	/*! \brief unpack a vector
	 *
	 * \warning the properties should match the packed properties, and the obj must have the same size of the packed vector, consider to pack
	 *          this information if you cannot infer-it
	 *
	 * \param ext preallocated memory from where to unpack the vector
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp> void unpack(ExtPreAlloc<HeapMemory> & mem, openfpm::vector<T> & obj, Unpack_stat & ps)
	{
		// if no properties should be unpacked return
		if (sizeof...(prp) == 0)
			return;

		size_t id = 0;

		// Sending property object
		typedef openfpm::vector<T> vctr;
		typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> stype;
		stype svect;

		// Calculate the size to pack the object
		size_t size = svect.calculateMem(obj.size(),0);

		// Create a Pointer object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());
		auto obj_it = obj.getIterator();

		while (obj_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_dst;
			// destination object type
			typedef encapc<1,prp_object,typename stype::memory_conf > encap_src;

			// Copy only the selected properties
			object_s_di<encap_src,encap_dst,ENCAP,prp...>(src.get(id),obj.get(obj_it.get()));

			++id;
			++obj_it;
		}

		ps.addOffset(size);
	}


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
	 * \param v element to add
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
#ifdef DEBUG
		if (key >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << key << "\n";
		}
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
	 * \tparam p must be 0
	 *
	 * \param id element to get
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline const T& get(size_t id) const
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
	 * \param id element to get
	 *
	 * \return the element reference
	 *
	 */
	inline T & get(size_t id)
	{
#ifdef DEBUG
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
#ifdef DEBUG
		if (id >= base.size())
		{
			std::cerr << "Error vector: " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";
		}
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
	vector() {}

	//! Constructor, vector of size sz
	vector(size_t sz):base(sz) {}

	/*! swap the content of the vector
	 *
	 * \param v vector to be swapped with
	 *
	 */
	void swap(openfpm::vector<T,HeapMemory,grow_policy_double,STD_VECTOR> & v)
	{
		base.swap(v.base);
	}

	/*! \brief Get iterator
	 *
	 * \return an iterator
	 *
	 */

	vector_key_iterator getIterator() const
	{
		return vector_key_iterator(base.size());
	}

	// Structures to check if object has or has not calculateMem() function

	template<bool cond, typename T1>
	struct calculateMem_cond
	{
		size_t calculateMemory(T1 & obj, size_t n, size_t e)
		{
			return grow_policy_double::grow(0,n) * sizeof(T1);
		}

	};

	template<typename T1>
	struct calculateMem_cond<true, T1>
	{
		size_t calculateMemory(T1 & obj, size_t n, size_t e)
		{
			size_t res = 0;
			size_t count = grow_policy_double::grow(0,n) - obj.size();
			for (int i = 0; i < n; i++) {
				res += obj.get(i).calculateMem(n,0);
			}
			return res+count*sizeof(T1);
		}
	};

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
	inline size_t calculateMem(size_t n, size_t e)
	{
		if (n == 0)
			return 0;

		calculateMem_cond<has_calculateMem<T>::type::value, openfpm::vector<T, HeapMemory, grow_policy_double>> cm;
		return cm.calculateMemory(*this,n,0);
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

	/*! \brief This class has pointer inside
	 *
	 * \return false
	 *
	 */
	static bool noPointers()
	{
		return false;
	}
};


#endif /* MAP_VECTOR_STD_HPP_ */
